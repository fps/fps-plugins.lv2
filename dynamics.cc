#include "generated/ttl2c_relative_dynamics.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

typedef struct plugin_state
{
    float abs1;
    float abs2;
    float samplerate;
    std::vector<float> buffer;
    int buffer_head;

} plugin_state_t;

static plugin_t* instantiate(plugin_t *instance, double sample_rate, const char *bundle_path, const LV2_Feature *const *features)
{
    instance->state = new plugin_state_t;
    instance->state->samplerate = sample_rate;
    instance->state->abs1 = 0;
    instance->state->abs2 = 0;
    instance->state->buffer = std::vector<float>(2*sample_rate, 0);
    instance->state->buffer_head = 0; 

    return instance;
}

static void cleanup(plugin_t *instance) {
    delete instance->state;
}

#define EPSILON 0.0001f

static void run(
    plugin_t *instance, uint32_t nframes, 
    const plugin_port_in_t in, const plugin_port_out_t out, 
    const plugin_port_t1_t t1, const plugin_port_t2_t t2, 
    const plugin_port_strength_t strength, const plugin_port_delay_t delay, 
    const plugin_port_maxratio_t maxratio, const plugin_port_minratio_t minratio
)
{
    plugin_state_t *tinstance = instance->state;

    const float a1 = 1.0f - expf((-1.0f/tinstance->samplerate) / (t1.data[0] / 1000.0f));
    const float a2 = 1.0f - expf((-1.0f/tinstance->samplerate) / (t2.data[0] / 1000.0f));

    for(uint32_t sample_index = 0; sample_index < nframes; ++sample_index)
    {
        tinstance->abs1 = a1 * fabs(in.data[sample_index]) + (1.0f - a1) * tinstance->abs1;
        tinstance->abs2 = a2 * fabs(in.data[sample_index]) + (1.0f - a2) * tinstance->abs2;

        const float r = (EPSILON + tinstance->abs1) / (EPSILON + tinstance->abs2);
        float scale = powf(1.0f / r, strength.data[0]);

        if (scale > maxratio.data[0]) {
            scale = maxratio.data[0];
        }

        if (scale < minratio.data[0]) {
            scale = minratio.data[0];
        }


        tinstance->buffer[tinstance->buffer_head] = in.data[sample_index];
        
        int buffer_tail = tinstance->buffer_head - tinstance->samplerate * (delay.data[0] / 1000);
        if (buffer_tail < 0) {
            buffer_tail += 2 * tinstance->samplerate;
        }
        
        out.data[sample_index] = scale * tinstance->buffer[buffer_tail];
        ++tinstance->buffer_head;
        tinstance->buffer_head %= 2 * (int)tinstance->samplerate;
    }
}

static const plugin_callbacks_t plugin_callbacks =
{
    .instantiate = instantiate,
    .run = run,
    .cleanup = cleanup,
};

#include "generated/ttl2c_relative_dynamics.c"

