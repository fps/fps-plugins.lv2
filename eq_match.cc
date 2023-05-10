#include "generated/ttl2c_eq_match.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <fftw3.h>
#include <iostream>
#include <complex>

#include <FFTConvolver/FFTConvolver.h>

#define FPS_PLUGINS_EQ_MATCH_FLOATING_POINT_TYPE float
#include "eq_match.h"

#define FFT_SIZE 2048

#define BLOCK_SIZE 32

typedef struct plugin_state
{
  plugin_state (float sample_rate, size_t fft_size) :
    m_buffer (100000),
    m_match (sample_rate, fft_size)
  {
    m_convolver.init (BLOCK_SIZE, m_match.m_minimum_phase_response, FFT_SIZE);
  }

  std::vector<float> m_buffer;

  eq_match m_match;

  bool m_previous_analyze1;
  bool m_previous_analyze2;

  fftconvolver::FFTConvolver m_convolver;

  std::vector<float> m_convolution_buffer;
  size_t m_convolution_buffer_head;
}
plugin_state_t;

static plugin_t* instantiate(plugin_t *instance, double sample_rate, const char *bundle_path, const LV2_Feature *const *features)
{
    instance->state = new plugin_state_t (FFT_SIZE, sample_rate);

    instance->state->m_convolution_buffer.resize(FFT_SIZE, 0);
    instance->state->m_convolution_buffer_head = FFT_SIZE - 1;

    instance->state->m_previous_analyze1 = false;
    instance->state->m_previous_analyze2 = false;

    return instance;
}

static void cleanup(plugin_t *instance)
{
    // fftw_cleanup ();
    delete instance->state;
}

#define EPSILON 0.0001f

static void run
(
    plugin_t *instance, uint32_t nframes, 
    const plugin_port_in_t in, const plugin_port_out_t out, 
    const plugin_port_analyze1_t analyze1, const plugin_port_analyze2_t analyze2, 
    const plugin_port_apply_t apply, const plugin_port_gain_t gain
)
{
    plugin_state_t *tinstance = instance->state;

    const bool previous_analyze1 = tinstance->m_previous_analyze1;
    const bool previous_analyze2 = tinstance->m_previous_analyze2;

    if (!previous_analyze1 && analyze1.data > 0)
    {
      tinstance->m_match.reset_buffer1 ();
    }

    if (!previous_analyze2 && analyze2.data > 0)
    {
      tinstance->m_match.reset_buffer2 ();
    }

    if ((previous_analyze1 && !(analyze1.data > 0)) || (previous_analyze2 && !(analyze2.data > 0)))
    {
      tinstance->m_convolver.reset ();
      tinstance->m_match.calculate_response ();
      tinstance->m_convolver.init (BLOCK_SIZE,  tinstance->m_match.m_minimum_phase_response, FFT_SIZE);
    }

    if (analyze1.data > 0)
    {
      tinstance->m_match.add_frames_to_buffer1 (in.data, nframes);
    }

    if (analyze2.data > 0)
    {
      tinstance->m_match.add_frames_to_buffer2 (in.data, nframes);
    }

    /*
    for(uint32_t sample_index = 0; sample_index < nframes; ++sample_index)
    {
      tinstance->m_convolution_buffer[tinstance->m_convolution_buffer_head] = in.data[sample_index];

      out.data[sample_index] = in.data[sample_index];

      if (apply.data > 0)
      {
          out.data[sample_index] = 0;
          for (size_t index = 0; index < FFT_SIZE; ++index)
          {
              out.data[sample_index] += tinstance->m_match.m_minimum_phase_response[index] * tinstance->m_convolution_buffer[(tinstance->m_convolution_buffer_head + index) % FFT_SIZE];
          }
      }

      if (tinstance->m_convolution_buffer_head == 0)
      {
        tinstance->m_convolution_buffer_head = FFT_SIZE - 1;
      }
      else
      {
        --tinstance->m_convolution_buffer_head;
      }
    }
    */

    if (apply.data > 0)
    {
      tinstance->m_convolver.process (in.data, &tinstance->m_buffer[0], nframes);
      memcpy (out.data, &tinstance->m_buffer[0], nframes * sizeof(float));
    }
    else
    {
      memcpy (out.data, in.data, nframes * sizeof (float));
    }

    tinstance->m_previous_analyze1 = analyze1.data > 0;
    tinstance->m_previous_analyze2 = analyze2.data > 0;
}

static const plugin_callbacks_t plugin_callbacks = {
    .instantiate = instantiate,
    .run = run,
    .cleanup = cleanup,
};

#include "generated/ttl2c_eq_match.c"

