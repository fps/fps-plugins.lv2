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
#define EQ_MATCH_STATE_SAMPLE_RATE EQ_MATCH_URI "#sample_rate"

#define EQ_MATCH_FFT_TIME 0.05
#define EQ_MATCH_BLOCK_SIZE 32

// #define EQ_MATCH_LOG(x) { std::cerr << x; }
#define EQ_MATCH_LOG(x) { }

struct plugin_state
{
    eq_match m_match;

    bool m_previous_analyze1;
    bool m_previous_analyze2;

    fftconvolver::FFTConvolver m_linear_phase_convolver;
    fftconvolver::FFTConvolver m_minimum_phase_convolver;

    plugin_state (size_t fft_size) :
        m_match (fft_size),
        m_previous_analyze1 (false),
        m_previous_analyze2 (false)
    {
        reset ();
    }

    void reset ()
    {
        m_match.reset ();
        init ();
    }

    void init ()
    {
        m_linear_phase_convolver.reset ();
        m_minimum_phase_convolver.reset ();

        m_linear_phase_convolver.init
            (EQ_MATCH_BLOCK_SIZE, &m_match.m_linear_phase_response[0], m_match.m_linear_phase_response.size ());

        m_minimum_phase_convolver.init
            (EQ_MATCH_BLOCK_SIZE, &m_match.m_minimum_phase_response[0], m_match.m_minimum_phase_response.size ());
    }
};

struct plugin
{
    std::atomic<bool> m_working;
    LV2_Worker_Schedule m_worker_schedule;
    LV2_URID_Map m_urid_map;
    float m_sample_rate;
    
    std::vector<float *> m_ports;

    enum WORKER_COMMAND
    {
        ADD_FRAMES_TO_BUFFER1 = 1,
        ADD_FRAMES_TO_BUFFER2,
        CALCULATE_RESPONSE,
        RESET_BUFFER1,
        RESET_BUFFER2
    };

    std::vector<float> m_worker_schedule_buffer;

    plugin_state m_plugin_state;

    plugin (float sample_rate, size_t fft_size) :
        m_working (false),
        m_sample_rate (sample_rate),
        m_ports (8, 0),
        // see comment above
        m_worker_schedule_buffer (fft_size + 2, 0),
        m_plugin_state (fft_size)
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
    EQ_MATCH_LOG("eq_match: instantiate\n")
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

    plugin *instance = new plugin (sample_rate, (size_t)(EQ_MATCH_FFT_TIME * sample_rate));
    instance->m_worker_schedule = worker_schedule;
    instance->m_urid_map = urid_map;
    return (LV2_Handle)instance;
}

static void activate (LV2_Handle instance)
{
    EQ_MATCH_LOG("eq_match: activate\n")
    ((plugin*)instance)->m_plugin_state.init ();
}

static void cleanup(LV2_Handle instance)
{
    EQ_MATCH_LOG("eq_match: cleanup\n")
    // fftw_cleanup ();
    delete ((plugin*)instance);
}

static void connect_port (LV2_Handle instance, uint32_t port, void *data_location)
{
        plugin &p = *((plugin*)instance);
        p.m_ports[port] = (float*)data_location;
}

static void schedule_worker_command (plugin &the_plugin, int c1, int c2 = 0, float *data = 0, size_t data_size = 0)
{
    assert (data_size <= the_plugin.m_plugin_state.m_match.m_fft_size);

    the_plugin.m_worker_schedule_buffer[0] = (float)c1;
    the_plugin.m_worker_schedule_buffer[1] = (float)c2;

    for (size_t index = 0; index < data_size; ++index)
    {
        the_plugin.m_worker_schedule_buffer[2 + index] = data[index];
    }

    LV2_Worker_Status status = the_plugin.m_worker_schedule.schedule_work
    (
        the_plugin.m_worker_schedule.handle,
        (2 + data_size) * sizeof (float),
        &the_plugin.m_worker_schedule_buffer[0]
    );
    
    if (status != LV2_WORKER_SUCCESS)
    {
        std::cerr << "Failed to schedule work: " << c1 << " " << c2 << "\n;";
    }
}

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
    const float &response_gain = *the_plugin.m_ports[6];
    const float &overall_gain  = *the_plugin.m_ports[7];

    const bool &previous_analyze1 = state.m_previous_analyze1;
    const bool &previous_analyze2 = state.m_previous_analyze2;

    if (!previous_analyze1 && analyze1 > 0)
    {
        schedule_worker_command (the_plugin, plugin::WORKER_COMMAND::RESET_BUFFER1);
    }

    if (!previous_analyze2 && analyze2 > 0)
    {
        schedule_worker_command (the_plugin, plugin::WORKER_COMMAND::RESET_BUFFER2);
    }

    if ((previous_analyze1 && !(analyze1 > 0)) || (previous_analyze2 && !(analyze2 > 0)))
    {
        the_plugin.m_working = true;
        schedule_worker_command (the_plugin, plugin::WORKER_COMMAND::CALCULATE_RESPONSE);
    }

    const size_t fft_size = state.m_match.m_fft_size;

    if (analyze1 > 0)
    {
        size_t samples_remaining = sample_count;
        size_t samples_written = 0;
        while (samples_remaining)
        {
            the_plugin.m_worker_schedule_buffer[0] =
                plugin::WORKER_COMMAND::ADD_FRAMES_TO_BUFFER1;

            size_t samples_in_buffer = 0;
            for
            (
                size_t index = 0;
                samples_remaining > 0 && index < fft_size;
                ++index
            )
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
            the_plugin.m_worker_schedule_buffer[0] =
                plugin::WORKER_COMMAND::ADD_FRAMES_TO_BUFFER2;

            size_t samples_in_buffer = 0;
            for
            (
                size_t index = 0;
                samples_remaining > 0 && index < fft_size;
                ++index
            )
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
        
        const float response_gain_factor = pow(10, response_gain/20);
        for (size_t index = 0; index < sample_count; ++index)
        {
            out[index] *= response_gain_factor;
        }
    }
    else
    {
        memcpy (out, in, sample_count * sizeof (float));
    }

    const float overall_gain_factor = pow(10, overall_gain/20);
    for (size_t index = 0; index < sample_count; ++index)
    {
        out[index] *= overall_gain_factor;
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
    EQ_MATCH_LOG("eq_match: save state\n")
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

    LV2_URID sample_rate_urid =
        the_plugin.m_urid_map.map
        (
            the_plugin.m_urid_map.handle,
            EQ_MATCH_STATE_SAMPLE_RATE
        );

    LV2_State_Status status;

    plugin_state &the_state = the_plugin.m_plugin_state;

    status = store
    (
        handle,
        linear_phase_response_urid,
        &the_state.m_match.m_linear_phase_response[0],
        the_state.m_match.m_linear_phase_response.size () * sizeof (float),
        the_plugin.m_urid_map.map (the_plugin.m_urid_map.handle, LV2_ATOM__Chunk),
        LV2_STATE_IS_POD
    );

    if (status != LV2_STATE_SUCCESS) 
    {
        std::cerr << "eq_match: Failed to save linear phase response\n";
        return LV2_STATE_ERR_UNKNOWN;
    }

    status = store
    (
        handle,
        minimum_phase_response_urid,
        &the_state.m_match.m_minimum_phase_response[0],
        the_state.m_match.m_minimum_phase_response.size () * sizeof (float),
        the_plugin.m_urid_map.map (the_plugin.m_urid_map.handle, LV2_ATOM__Chunk),
        LV2_STATE_IS_POD
    );

    if (status != LV2_STATE_SUCCESS)
    {
        std::cerr << "eq_match: Failed to save minimum phase response\n";
        return LV2_STATE_ERR_UNKNOWN;
    }

    status = store
    (
        handle,
        sample_rate_urid,
        &the_plugin.m_sample_rate,
        sizeof (float),
        the_plugin.m_urid_map.map (the_plugin.m_urid_map.handle, LV2_ATOM__Float),
        LV2_STATE_IS_POD
    );
        
    if (status != LV2_STATE_SUCCESS)
    {
        std::cerr << "eq_match: Failed to save sample rate\n";
        return LV2_STATE_ERR_UNKNOWN;
    }

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
    EQ_MATCH_LOG("eq_match: restore state\n")

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
    uint32_t type;
    uint32_t the_flags;

    EQ_MATCH_LOG("retrieving state 1...\n")

    plugin_state &the_state = the_plugin.m_plugin_state;

    float *linear_phase_response_data =
        (float*)retrieve
        (
            handle,
            linear_phase_response_urid,
            &size,
            &type,
            &the_flags
        );

    if (size != sizeof (float) * the_state.m_match.m_linear_phase_response.size ()) return LV2_STATE_ERR_UNKNOWN;
    if (linear_phase_response_data == 0) return LV2_STATE_ERR_UNKNOWN;

    EQ_MATCH_LOG("retrieving state 2...\n")

    float *minimum_phase_response_data =
        (float*)retrieve
        (
            handle,
            minimum_phase_response_urid,
            &size,
            &type,
            &the_flags
        );

    if (size != sizeof (float) * the_state.m_match.m_minimum_phase_response.size ()) return LV2_STATE_ERR_UNKNOWN;
    if (minimum_phase_response_data == 0) return LV2_STATE_ERR_UNKNOWN;

    EQ_MATCH_LOG("setting up convolvers...\n")

    the_state.m_match.reset ();
    std::copy
    (
        linear_phase_response_data, 
        linear_phase_response_data + the_state.m_match.m_linear_phase_response.size (), 
        the_state.m_match.m_linear_phase_response.begin ()
    );
    
    std::copy
    (
        minimum_phase_response_data, 
        minimum_phase_response_data + the_state.m_match.m_minimum_phase_response.size (), 
        the_state.m_match.m_minimum_phase_response.begin ()
    );

    the_state.init ();
    
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
    EQ_MATCH_LOG("eq_match: work: " << size << "\n")

    plugin &the_plugin = *((plugin*)instance);
    plugin_state &the_state = the_plugin.m_plugin_state;
    eq_match &the_match = the_state.m_match;

    if (size >= 2 * sizeof (float))
    {
        float the_cmd = ((float*)data)[0];
        float the_size = ((float*)data)[1];
        EQ_MATCH_LOG(the_cmd << " " << the_size << "\n")

        if (the_cmd == (float)plugin::WORKER_COMMAND::ADD_FRAMES_TO_BUFFER1)
        {
            EQ_MATCH_LOG("work: add frames to buffer 1: " << the_size << " frames.\n")
            the_match.add_frames_to_buffer1 (((float*)data)+2, the_size);
        }
        else
        if (the_cmd == (float)plugin::WORKER_COMMAND::ADD_FRAMES_TO_BUFFER2)
        {
            EQ_MATCH_LOG("work: add frames to buffer 2: " << the_size << " frames.\n")
            the_match.add_frames_to_buffer2 (((float*)data)+2, the_size);
        }
        else
        if (the_cmd == (float)plugin::WORKER_COMMAND::CALCULATE_RESPONSE)
        {
            EQ_MATCH_LOG("work: calculate response\n")
            the_match.calculate_response ();

            the_state.init ();

            the_plugin.m_working = false;
            EQ_MATCH_LOG("work: done.\n")
        }
        else
        if (the_cmd == (float)plugin::WORKER_COMMAND::RESET_BUFFER1)
        {
            the_match.reset_buffer1 ();
        }
        else
        if (the_cmd == (float)plugin::WORKER_COMMAND::RESET_BUFFER2)
        {
            the_match.reset_buffer2 ();
        }
        else
        {
            EQ_MATCH_LOG("eq_match: Unknown command\n")
            return LV2_WORKER_ERR_UNKNOWN;
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
    EQ_MATCH_LOG("get extension_data. URI: " << uri << "\n")
    if (std::string(uri) == LV2_STATE__interface)
    {
        EQ_MATCH_LOG("get extension data for state\n")
        return &state_interface;
    }

    if (std::string(uri) == LV2_WORKER__interface)
    {
        EQ_MATCH_LOG("get extension data for worker\n")
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

