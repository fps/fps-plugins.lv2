#include <atomic>

#include <lv2.h>
#include <lv2/worker/worker.h>
#include <lv2/state/state.h>
#include <lv2/urid/urid.h>
#include <lv2/atom/atom.h>

#include <FFTConvolver/FFTConvolver.h>

#include "eq_match.h"

#include "common.h"

#define EQ_MATCH_URI FPS_PLUGINS_BASE_URI "/eq_match"
#define EQ_MATCH_STATE_URI EQ_MATCH_URI "#state"
#define EQ_MATCH_STATE_LINEAR_PHASE_RESPONSE EQ_MATCH_URI "#linear_phase_response"
#define EQ_MATCH_STATE_MINIMUM_PHASE_RESPONSE EQ_MATCH_URI "#minimal_phase_response"

#define FFT_SIZE 2048
#define BLOCK_SIZE 32

struct plugin_state
{
    eq_match m_match;

    bool m_previous_analyze1;
    bool m_previous_analyze2;

    fftconvolver::FFTConvolver m_linear_phase_convolver;
    fftconvolver::FFTConvolver m_minimum_phase_convolver;

    plugin_state (float sample_rate, size_t fft_size) :
        m_match (sample_rate, fft_size),
        m_previous_analyze1 (false),
        m_previous_analyze2 (false)
    {
      reset ();
    }

    void reset ()
    {
      m_match.reset ();

      m_linear_phase_convolver.init
          (BLOCK_SIZE, &m_match.m_minimum_phase_response[0], FFT_SIZE);

      m_minimum_phase_convolver.init
          (BLOCK_SIZE, &m_match.m_minimum_phase_response[0], FFT_SIZE);
    }
};

struct plugin
{
    size_t m_fft_size;
    std::atomic<bool> m_working;
    LV2_Worker_Schedule m_worker_schedule;
    LV2_URID_Map m_urid_map;

    std::vector<float *> m_ports;

    // first entry is a command. second entry is the number of samples,
    // then come the samples. Commands:
    // 1 - add frames to buffer 1
    // 2 - add frames to buffer 2
    // 3 - calculate calculate (second entry == 0)
    // 4 - reset (second entry == 0)
    // 5 - linear phase response (work response)
    // 6 - minimal phase response (work response)
    std::vector<float> m_worker_schedule_buffer;

    plugin_state m_plugin_state;

    plugin (float sample_rate, size_t fft_size) :
        m_fft_size (fft_size),
        m_working (false),
        m_ports (7, 0),
        // see comment above
        m_worker_schedule_buffer (fft_size + 2, 0),
        m_plugin_state (sample_rate, fft_size)
    {

    }
};

LV2_Handle instantiate
(
    const LV2_Descriptor *descriptor,
    double sample_rate,
    const char *bundle_path,
    const LV2_Feature *const *features
)
{
    const LV2_Feature * feature = 0;
    bool worker_found = false;
    LV2_Worker_Schedule worker_schedule;
    LV2_URID_Map urid_map;

    bool urid_map_found = false;
    while ((feature = *(++features)) != 0)
    {
        if (std::string (feature->URI) == LV2_WORKER__schedule)
        {
            worker_found = true;
            worker_schedule = *((LV2_Worker_Schedule*)feature->data);
        }

        if (std::string (feature->URI) == LV2_URID__map)
        {
            urid_map_found = true;
            urid_map = *((LV2_URID_Map*)feature->data);
        }
    }

    if (false == worker_found) return 0;
    if (false == urid_map_found) return 0;

    plugin *instance = new plugin (FFT_SIZE, sample_rate);
    instance->m_worker_schedule = worker_schedule;
    instance->m_urid_map = urid_map;
    return (LV2_Handle)instance;
}

static void activate (LV2_Handle instance)
{
  ((plugin*)instance)->m_plugin_state.reset ();
}

static void cleanup(LV2_Handle instance)
{
    // fftw_cleanup ();
    delete ((plugin*)instance);
}

static void connect_port (LV2_Handle instance, uint32_t port, void *data_location)
{
        plugin &p = *((plugin*)instance);
        p.m_ports[port] = (float*)data_location;
}

// #define EPSILON 0.0001f

static void run
(
    LV2_Handle instance,
    uint32_t sample_count
)
{
    plugin &the_plugin = *((plugin*)instance);
    plugin_state &state = the_plugin.m_plugin_state;

    // Audio ports
    const float *in            =  the_plugin.m_ports[0];
    float       *out           =  the_plugin.m_ports[1];

    // Control eports
    const float &analyze1      = *the_plugin.m_ports[2];
    const float &analyze2      = *the_plugin.m_ports[3];
    const float &apply         = *the_plugin.m_ports[4];
    const float &minimum_phase = *the_plugin.m_ports[5];
    const float &gain          = *the_plugin.m_ports[6];

    const bool &previous_analyze1 = state.m_previous_analyze1;
    const bool &previous_analyze2 = state.m_previous_analyze2;

    if (!previous_analyze1 && analyze1 > 0)
    {
        state.m_match.reset_buffer1 ();
    }

    if (!previous_analyze2 && analyze2 > 0)
    {
        state.m_match.reset_buffer2 ();
    }

    if ((previous_analyze1 && !(analyze1 > 0)) || (previous_analyze2 && !(analyze2 > 0)))
    {
        the_plugin.m_working = true;

        the_plugin.m_worker_schedule_buffer[0] = 3.0f;
        the_plugin.m_worker_schedule_buffer[1] = 0;
        the_plugin.m_worker_schedule.schedule_work
        (
            the_plugin.m_worker_schedule.handle,
            2 * sizeof (float),
            &the_plugin.m_worker_schedule_buffer[0]
        );
    }

    if (analyze1 > 0)
    {
        size_t samples_remaining = sample_count;
        size_t samples_written = 0;
        while (samples_remaining)
        {
            the_plugin.m_worker_schedule_buffer[0] = 1.0f;
            size_t samples_in_buffer = 0;
            for (size_t index = 0; samples_remaining > 0 && index < the_plugin.m_fft_size; ++index)
            {
                the_plugin.m_worker_schedule_buffer[2+index] = in[samples_written];
                ++samples_in_buffer;
                ++samples_written;
                --samples_remaining;
            }

            the_plugin.m_worker_schedule_buffer[1] = samples_in_buffer;

            the_plugin.m_worker_schedule.schedule_work
            (
                the_plugin.m_worker_schedule.handle,
                (2 + samples_in_buffer) * sizeof (float),
                &the_plugin.m_worker_schedule_buffer[0]
            );
        }
    }

    if (analyze2 > 0)
    {
        size_t samples_remaining = sample_count;
        size_t samples_written = 0;
        while (samples_remaining)
        {
            the_plugin.m_worker_schedule_buffer[0] = 2.0f;
            size_t samples_in_buffer = 0;
            for (size_t index = 0; samples_remaining > 0 && index < the_plugin.m_fft_size; ++index)
            {
                the_plugin.m_worker_schedule_buffer[2+index] = in[samples_written];
                ++samples_in_buffer;
                ++samples_written;
                --samples_remaining;
            }

            the_plugin.m_worker_schedule_buffer[1] = samples_in_buffer;

            the_plugin.m_worker_schedule.schedule_work
            (
                the_plugin.m_worker_schedule.handle,
                (2 + samples_in_buffer) * sizeof (float),
                &the_plugin.m_worker_schedule_buffer[0]
            );
        }
    }

    if (apply > 0 && !the_plugin.m_working)
    {
        if (minimum_phase > 0)
        {
            state.m_minimum_phase_convolver.process (in, out, sample_count);
        }
        else
        {
            state.m_linear_phase_convolver.process (in, out, sample_count);
        }
    }
    else
    {
        memcpy (out, in, sample_count * sizeof (float));
    }

    const float gain_factor = pow(10, gain/20);
    for (size_t index = 0; index < sample_count; ++index)
    {
        out[index] *= gain_factor;
    }

    state.m_previous_analyze1 = analyze1 > 0;
    state.m_previous_analyze2 = analyze2 > 0;
}

static LV2_State_Status save_state
(
    LV2_Handle instance,
    LV2_State_Store_Function store,
    LV2_State_Handle handle,
    uint32_t flags,
    const LV2_Feature *const *features
)
{
    std::cerr << "save state\n";
    // plugin_t *pinstance = ((plugin_t*)instance);
    // std::cerr << pinstance->map (
    //store (handle,
    plugin &the_plugin = *((plugin*)instance);

    LV2_URID linear_phase_response_urid =
        the_plugin.m_urid_map.map
        (
            the_plugin.m_urid_map.handle,
            EQ_MATCH_STATE_LINEAR_PHASE_RESPONSE
        );

    LV2_URID minimum_phase_response_urid =
        the_plugin.m_urid_map.map
        (
            the_plugin.m_urid_map.handle,
            EQ_MATCH_STATE_MINIMUM_PHASE_RESPONSE
        );

    LV2_State_Status status;

    status = store
    (
        handle,
        linear_phase_response_urid,
        &the_plugin.m_plugin_state.m_match.m_linear_phase_response[0],
        the_plugin.m_fft_size * sizeof (float),
        the_plugin.m_urid_map.map (the_plugin.m_urid_map.handle, LV2_ATOM__Chunk),
        LV2_STATE_IS_POD
    );

    if (status != LV2_STATE_SUCCESS) return LV2_STATE_ERR_UNKNOWN;

    status = store
    (
        handle,
        minimum_phase_response_urid,
        &the_plugin.m_plugin_state.m_match.m_minimum_phase_response[0],
        the_plugin.m_fft_size * sizeof (float),
        the_plugin.m_urid_map.map (the_plugin.m_urid_map.handle, LV2_ATOM__Chunk),
        LV2_STATE_IS_POD
    );

    if (status != LV2_STATE_SUCCESS) return LV2_STATE_ERR_UNKNOWN;

    return LV2_STATE_SUCCESS;
}

static LV2_State_Status restore_state
(
    LV2_Handle instance,
    LV2_State_Retrieve_Function retrieve,
    LV2_State_Handle handle,
    uint32_t flags,
    const LV2_Feature *const *features
)
{
    std::cerr << "restore state\n";

    plugin &the_plugin = *((plugin*)instance);

    LV2_URID linear_phase_response_urid =
        the_plugin.m_urid_map.map
        (
            the_plugin.m_urid_map.handle,
            EQ_MATCH_STATE_LINEAR_PHASE_RESPONSE
        );

    LV2_URID minimum_phase_response_urid =
        the_plugin.m_urid_map.map
        (
            the_plugin.m_urid_map.handle,
            EQ_MATCH_STATE_MINIMUM_PHASE_RESPONSE
        );

    size_t size;

    std::cerr << "retrieving state 1...\n";

    float *linear_phase_response_data =
        (float*)retrieve
        (
            handle,
            linear_phase_response_urid,
            &size,
            0,
            0
        );

    if (size != sizeof (float) * the_plugin.m_fft_size) return LV2_STATE_ERR_UNKNOWN;

    std::cerr << "retrieving state 2...\n";

    float *minimum_phase_response_data =
        (float*)retrieve
        (
            handle,
            minimum_phase_response_urid,
            &size,
            0,
            0
        );

    if (size != sizeof (float) * the_plugin.m_fft_size) return LV2_STATE_ERR_UNKNOWN;

    std::cerr << "setting up convolvers...\n";

    the_plugin.m_plugin_state.m_match.reset ();

    the_plugin.m_plugin_state.m_linear_phase_convolver.init
        (BLOCK_SIZE, linear_phase_response_data, FFT_SIZE);

    the_plugin.m_plugin_state.m_minimum_phase_convolver.init
        (BLOCK_SIZE, minimum_phase_response_data, FFT_SIZE);

    return LV2_STATE_SUCCESS;
}

static LV2_State_Interface state_interface =
{
    .save = save_state,
    .restore = restore_state,
};

LV2_Worker_Status work
(
    LV2_Handle instance,
    LV2_Worker_Respond_Function respond,
    LV2_Worker_Respond_Handle handle,
    uint32_t size,
    const void *data
)
{
    plugin &the_plugin = *((plugin*)instance);

    // std::cerr << "work: " << size << "\n";
    if (size >= 2 * sizeof (float))
    {
        float the_cmd = ((float*)data)[0];
        float the_size = ((float*)data)[1];
        // std::cerr << the_cmd << " " << the_size << "\n";

        if (the_cmd == 1)
        {
            // std::cerr << "work: add frames to buffer 1: " << the_size << " frames.\n";
            the_plugin.m_plugin_state.m_match.add_frames_to_buffer1 (((float*)data)+2, the_size);
        }
        else
        if (the_cmd == 2)
        {
            // std::cerr << "work: add frames to buffer 2: " << the_size << " frames.\n";
            the_plugin.m_plugin_state.m_match.add_frames_to_buffer2 (((float*)data)+2, the_size);
        }
        else
        if (the_cmd == 3)
        {
            // std::cerr << "work: calculate response\n";
            the_plugin.m_plugin_state.m_linear_phase_convolver.reset ();
            the_plugin.m_plugin_state.m_minimum_phase_convolver.reset ();

            the_plugin.m_plugin_state.m_match.calculate_response ();

            the_plugin.m_plugin_state.m_linear_phase_convolver.init
                (BLOCK_SIZE, &the_plugin.m_plugin_state.m_match.m_linear_phase_response[0], FFT_SIZE);

            the_plugin.m_plugin_state.m_minimum_phase_convolver.init
                (BLOCK_SIZE, &the_plugin.m_plugin_state.m_match.m_minimum_phase_response[0], FFT_SIZE);

            the_plugin.m_working = false;
            // std::cerr << "work: done.\n";
        }
    }
    else
    {
        return LV2_WORKER_ERR_UNKNOWN;
    }
    return LV2_WORKER_SUCCESS;
}

LV2_Worker_Status work_response
(
    LV2_Handle instance,
    uint32_t size,
    const void *body
)
{
    return LV2_WORKER_SUCCESS;
}

LV2_Worker_Status end_run(LV2_Handle instance)
{
    return LV2_WORKER_SUCCESS;
}

static LV2_Worker_Interface worker_interface =
{
    .work = work,
    .work_response = work_response,
    .end_run = end_run,
};

static const void *extension_data (const char *uri)
{
    // std::cerr << "get extension_data. URI: " << uri << "\n";
    if (std::string(uri) == LV2_STATE__interface)
    {
        // std::cerr << "get extension data for state\n";
        return &state_interface;
    }

    if (std::string(uri) == LV2_WORKER__interface)
    {
        // std::cerr << "get extension data for worker\n";
        return &worker_interface;
    }

    return 0;
}

static LV2_Descriptor plugin_descriptor = {
    EQ_MATCH_URI,
    instantiate,
    connect_port,
    activate,
    run,
    0,
    cleanup,
    extension_data
};

LV2_SYMBOL_EXPORT const LV2_Descriptor* lv2_descriptor (uint32_t index) {
    if (0 == index) {
          return &plugin_descriptor;
    } else {
        return NULL;
    }
}

