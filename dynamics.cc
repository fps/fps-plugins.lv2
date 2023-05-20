#include <cmath>
#include <vector>
#include <lv2.h>

#include "common.h"

#define RELATIVE_DYNAMICS_URI FPS_PLUGINS_BASE_URI "/relative_dynamics"

struct plugin_state
{
    float m_sample_rate;
    float m_abs1;
    float m_abs2;
    std::vector<float> m_buffer;
    int m_buffer_head;

    plugin_state (float sample_rate) :
        m_sample_rate (sample_rate),
        m_buffer (2 * sample_rate, 0)
    {
        reset ();
    }

    void reset ()
    {
        m_abs1 = 0;
        m_abs2 = 0;
        std::fill (m_buffer.begin (), m_buffer.end (), 0.f);
        m_buffer_head = 0;
    }
};

struct plugin
{
    std::vector<float *> m_ports;

    plugin_state m_plugin_state;

    plugin (float sample_rate) :
        m_ports (8, 0),
        m_plugin_state (sample_rate)
    {

    }
};

static void connect_port (LV2_Handle instance, uint32_t port, void *data_location)
{
    plugin &p = *((plugin*)instance);
    p.m_ports[port] = (float*)data_location;
}

LV2_Handle instantiate
(
    const struct LV2_Descriptor *descriptor,
    double sample_rate,
    const char *bundle_path,
    const LV2_Feature *const *features
)
{
    plugin *instance = new plugin (sample_rate);
    return (LV2_Handle)instance;
}

static void activate (LV2_Handle instance)
{
    ((plugin*)instance)->m_plugin_state.reset ();
}

static void cleanup(LV2_Handle instance)
{
    delete (plugin*)instance;
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
    const float *in       =  the_plugin.m_ports[0];
    float       *out      =  the_plugin.m_ports[1];

    // Control eports
    const float &t1       = *the_plugin.m_ports[2];
    const float &t2       = *the_plugin.m_ports[3];
    const float &strength = *the_plugin.m_ports[4];
    const float &delay    = *the_plugin.m_ports[5];
    const float &maxratio = *the_plugin.m_ports[6];
    const float &minratio = *the_plugin.m_ports[7];

    const float a1 = 1.0f - expf((-1.0f/state.m_sample_rate) / (t1 / 1000.0f));
    const float a2 = 1.0f - expf((-1.0f/state.m_sample_rate) / (t2 / 1000.0f));

    for(uint32_t sample_index = 0; sample_index < sample_count; ++sample_index)
    {
        state.m_abs1 = a1 * fabs(in[sample_index]) + (1.0f - a1) * state.m_abs1;
        state.m_abs2 = a2 * fabs(in[sample_index]) + (1.0f - a2) * state.m_abs2;

        const float r = (EPSILON + state.m_abs1) / (EPSILON + state.m_abs2);
        float scale = powf(1.0f / r, strength);

        scale = std::max (std::min (scale, maxratio), minratio);

        state.m_buffer[state.m_buffer_head] = in[sample_index];
        
        int buffer_tail = state.m_buffer_head - state.m_sample_rate * (delay / 1000);
        if (buffer_tail < 0) 
        {
            buffer_tail += 2 * state.m_sample_rate;
        }
        
        out[sample_index] = scale * state.m_buffer[buffer_tail];
        ++state.m_buffer_head;
        state.m_buffer_head %= state.m_buffer.size ();
    }
}

static LV2_Descriptor plugin_descriptor = {
    RELATIVE_DYNAMICS_URI,
    instantiate,
    connect_port,
    activate,
    run,
    0,
    cleanup,
    0
};

LV2_SYMBOL_EXPORT const LV2_Descriptor* lv2_descriptor (uint32_t index) {
    if (0 == index) {
        return &plugin_descriptor;
    } else {
        return NULL;
    }
}

