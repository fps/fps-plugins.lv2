#include "generated/ttl2c_eq_match.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <fftw3.h>

#define FFT_SIZE 1024

typedef struct plugin_state {
    float samplerate;
    std::vector<float> hann_window;
    std::vector<float> buffer1;
    int buffer_head1;
    std::vector<float> buffer2;
    int buffer_head2;
    std::vector<float> spectrum1;
    std::vector<float> spectrum2;
    std::vector<float> response;
    bool previous_analyze1;
    bool previous_analyze2;

} plugin_state_t;

static plugin_t* instantiate(plugin_t *instance, double sample_rate, const char *bundle_path, const LV2_Feature *const *features) {
    instance->state = new plugin_state_t;
    instance->state->samplerate = sample_rate;
    instance->state->hann_window = std::vector<float>(FFT_SIZE, 0);
    instance->state->buffer1 = std::vector<float>(FFT_SIZE, 0);
    instance->state->buffer_head2 = 0;
    instance->state->buffer2 = std::vector<float>(FFT_SIZE, 0);
    instance->state->buffer_head2 = 0; 
    instance->state->spectrum1 = std::vector<float>(FFT_SIZE, 0);
    instance->state->spectrum2 = std::vector<float>(FFT_SIZE, 0);
    instance->state->response = std::vector<float>(FFT_SIZE, 0);
    instance->state->previous_analyze1 = false;
    instance->state->previous_analyze2 = false;

    for (size_t index = 0; index < FFT_SIZE; ++index) {
      instance->state->hann_window[index] = powf(sin(M_PI * index / FFT_SIZE), 2);
    }

    return instance;
}

static void cleanup(plugin_t *instance) {
    delete instance->state;
}

#define EPSILON 0.0001f

static void run(
    plugin_t *instance, uint32_t nframes, 
    const plugin_port_in_t in, const plugin_port_out_t out, 
    const plugin_port_analyze1_t analyze1, const plugin_port_analyze2_t analyze2, 
    const plugin_port_apply_t apply) 
{
    plugin_state_t *tinstance = instance->state;

    const bool previous_analyze1 = tinstance->previous_analyze1;
    const bool previous_analyze2 = tinstance->previous_analyze2;

    if (!previous_analyze1 && analyze1.data) {
      // start collecting samples for spectrum 1
      memset(&tinstance->buffer1[0], 0, FFT_SIZE);
    }

    if (!previous_analyze2 && analyze2.data) {
      // start collecting samples for spectrum 2
      memset(&tinstance->buffer2[0], 0, FFT_SIZE);
    }

    if ((previous_analyze1 && !analyze1.data) || (previous_analyze2 && !analyze2.data)) {
      // calculate response
    }

    for(uint32_t sample_index = 0; sample_index < nframes; ++sample_index) {
        tinstance->buffer1[tinstance->buffer_head1] = in.data[sample_index];
        
        int buffer_tail = tinstance->buffer_head1 - tinstance->samplerate * (delay.data / 1000);
        if (buffer_tail < 0) {
            buffer_tail += 2 * tinstance->samplerate;
        }
        
        out.data[sample_index] = scale * tinstance->buffer1[buffer_tail];
        ++tinstance->buffer_head1;
        tinstance->buffer_head1 %= 2 * (int)tinstance->samplerate;
    }

    tinstance->state->previous_analyze1 = analyze1;
    tinstance->state->previous_analyze2 = analyze2;
}

static const plugin_callbacks_t plugin_callbacks = {
    .instantiate = instantiate,
    .run = run,
    .cleanup = cleanup,
};

#include "generated/ttl2c_eq_match.c"

