#include <lv2.h>
#include <lv2/worker/worker.h>
#include <lv2/state/state.h>

#include <FFTConvolver/FFTConvolver.h>

#include "eq_match.h"

#include "common.h"
#define EQ_MATCH_URI FPS_PLUGINS_BASE_URI "/eq_match"
#define EQ_MATCH_STATE_URI EQ_MATCH_URI "#state"

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
          (BLOCK_SIZE, m_match.m_minimum_phase_response, FFT_SIZE);

      m_minimum_phase_convolver.init
          (BLOCK_SIZE, m_match.m_minimum_phase_response, FFT_SIZE);
    }
};

struct plugin
{
    std::vector<float *> m_ports;
    plugin_state m_plugin_state;

    plugin (float sample_rate, size_t fft_size) :
        m_ports (7, 0),
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
    plugin *instance = new plugin (FFT_SIZE, sample_rate);
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

#define EPSILON 0.0001f

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
        state.m_linear_phase_convolver.reset ();
        state.m_minimum_phase_convolver.reset ();

        state.m_match.calculate_response ();

        state.m_linear_phase_convolver.init
            (BLOCK_SIZE,    state.m_match.m_linear_phase_response, FFT_SIZE);

        state.m_minimum_phase_convolver.init
            (BLOCK_SIZE,    state.m_match.m_minimum_phase_response, FFT_SIZE);
    }

    if (analyze1 > 0)
    {
        state.m_match.add_frames_to_buffer1 (in, sample_count);
    }

    if (analyze2 > 0)
    {
        state.m_match.add_frames_to_buffer2 (in, sample_count);
    }

    if (apply > 0)
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

#define EQ_MATCH_STATE_SAMPLERATE "https://dfdx.eu/fps-plugins.lv2/eq_match#samplerate"

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
    std::cerr << "get extension_data. URI: " << uri << "\n";
    if (std::string(uri) == LV2_STATE__interface)
    {
        std::cerr << "get extension data for state\n";
        return &state_interface;
    }

    if (std::string(uri) == LV2_WORKER__interface)
    {
        std::cerr << "get extension data for worker\n";
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

