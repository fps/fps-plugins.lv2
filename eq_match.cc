#include "generated/ttl2c_eq_match.h"

#include <FFTConvolver/FFTConvolver.h>

#define EQ_MATCH_FLOAT float
#include "eq_match.h"

#define FFT_SIZE 2048

#define BLOCK_SIZE 32

typedef struct plugin_state
{
  plugin_state (float sample_rate, size_t fft_size) :
    m_match (sample_rate, fft_size),
    m_previous_analyze1 (false),
    m_previous_analyze2 (false)
  {
    m_linear_phase_convolver.init
      (BLOCK_SIZE, m_match.m_minimum_phase_response, FFT_SIZE);

    m_minimum_phase_convolver.init
      (BLOCK_SIZE, m_match.m_minimum_phase_response, FFT_SIZE);
  }

  eq_match m_match;

  bool m_previous_analyze1;
  bool m_previous_analyze2;

  fftconvolver::FFTConvolver m_linear_phase_convolver;
  fftconvolver::FFTConvolver m_minimum_phase_convolver;
}
plugin_state_t;

static plugin_t* instantiate
(
  plugin_t *instance, double sample_rate,
  const char *bundle_path,
  const LV2_Feature *const *features
)
{
    instance->state = new plugin_state_t (FFT_SIZE, sample_rate);
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
    const plugin_port_apply_t apply, const plugin_port_minimum_phase_t minimum_phase, const plugin_port_gain_t gain
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
      tinstance->m_linear_phase_convolver.reset ();
      tinstance->m_minimum_phase_convolver.reset ();

      tinstance->m_match.calculate_response ();

      tinstance->m_linear_phase_convolver.init
        (BLOCK_SIZE,  tinstance->m_match.m_linear_phase_response, FFT_SIZE);

      tinstance->m_minimum_phase_convolver.init
        (BLOCK_SIZE,  tinstance->m_match.m_minimum_phase_response, FFT_SIZE);
    }

    if (analyze1.data > 0)
    {
      tinstance->m_match.add_frames_to_buffer1 (in.data, nframes);
    }

    if (analyze2.data > 0)
    {
      tinstance->m_match.add_frames_to_buffer2 (in.data, nframes);
    }

    if (apply.data > 0)
    {
      if (minimum_phase.data > 0)
      {
        tinstance->m_minimum_phase_convolver.process (in.data, out.data, nframes);
      }
      else
      {
        tinstance->m_linear_phase_convolver.process (in.data, out.data, nframes);
      }
    }
    else
    {
      memcpy (out.data, in.data, nframes * sizeof (float));
    }

    const float gain_factor = pow(10, gain.data/20);
    for (size_t index = 0; index < nframes; ++index)
    {
      out.data[index] *= gain_factor;
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

