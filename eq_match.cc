#include "generated/ttl2c_eq_match.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <fftw3.h>
#include <iostream>

#define FFT_REAL_SIZE 1024
#define FFT_COMPLEX_SIZE FFT_REAL_SIZE

typedef struct plugin_state {
    float samplerate;
    std::vector<double> window;
    
    fftw_complex *buffer11;
    int buffer_head11;
    
    fftw_complex *buffer12;
    int buffer_head12;
    
    fftw_complex *buffer21;
    int buffer_head21;
    
    fftw_complex *buffer22;
    int buffer_head22;
    
    fftw_complex *spectrum1;
    fftw_complex *spectrum2;
    
    fftw_complex *response;
    
    bool previous_analyze1;
    bool previous_analyze2;
    
    fftw_plan fft_plan11;
    fftw_plan fft_plan12;
    fftw_plan fft_plan21;
    fftw_plan fft_plan22;
    
    fftw_plan ifft_plan;
    
    fftw_complex *fft_buffer1;
    size_t number_of_ffts1;
    
    fftw_complex *fft_buffer2;
    size_t number_of_ffts2;
    
    fftw_complex *fft_buffer3;
    fftw_complex *fft_buffer4;

    std::vector<float> convolution_buffer;
    size_t convolution_buffer_head;
} plugin_state_t;

static plugin_t* instantiate(plugin_t *instance, double sample_rate, const char *bundle_path, const LV2_Feature *const *features) {
    instance->state = new plugin_state_t;
    instance->state->samplerate = sample_rate;

    instance->state->window.resize (FFT_REAL_SIZE, 0);

    instance->state->buffer11 = (fftw_complex*)fftw_alloc_complex (FFT_COMPLEX_SIZE);
    instance->state->buffer_head11 = 0;
    instance->state->buffer12 = (fftw_complex*)fftw_alloc_complex (FFT_COMPLEX_SIZE);
    instance->state->buffer_head12 = FFT_REAL_SIZE/2; 
    
    instance->state->buffer21 = (fftw_complex*)fftw_alloc_complex (FFT_COMPLEX_SIZE);
    instance->state->buffer_head21 = 0;
    instance->state->buffer22 = (fftw_complex*)fftw_alloc_complex (FFT_COMPLEX_SIZE);
    instance->state->buffer_head22 = FFT_REAL_SIZE/2; 

    instance->state->spectrum1 = (fftw_complex*)fftw_alloc_complex (FFT_COMPLEX_SIZE);
    instance->state->spectrum2 = (fftw_complex*)fftw_alloc_complex (FFT_COMPLEX_SIZE);

    instance->state->response = (fftw_complex*)fftw_alloc_complex (FFT_COMPLEX_SIZE);

    instance->state->previous_analyze1 = false;
    instance->state->previous_analyze2 = false;

    instance->state->fft_buffer1 = (fftw_complex*)fftw_alloc_complex (FFT_COMPLEX_SIZE);
    instance->state->number_of_ffts1 = 0;
    
    instance->state->fft_buffer2 = (fftw_complex*)fftw_alloc_complex (FFT_COMPLEX_SIZE);
    instance->state->number_of_ffts2 = 0;
    
    instance->state->fft_buffer3 = (fftw_complex*)fftw_alloc_complex (FFT_COMPLEX_SIZE);
    instance->state->fft_buffer4 = (fftw_complex*)fftw_alloc_complex (FFT_COMPLEX_SIZE);

    for (size_t index = 0; index < FFT_REAL_SIZE; ++index) {
      instance->state->window[index] = pow(sin(M_PI * index / FFT_REAL_SIZE), 2);
    }

    instance->state->fft_plan11 = fftw_plan_dft_1d(FFT_REAL_SIZE, &instance->state->buffer11[0], instance->state->fft_buffer1, FFTW_FORWARD, FFTW_ESTIMATE);
    instance->state->fft_plan12 = fftw_plan_dft_1d(FFT_REAL_SIZE, &instance->state->buffer12[0], instance->state->fft_buffer1, FFTW_FORWARD, FFTW_ESTIMATE);
    
    instance->state->fft_plan21 = fftw_plan_dft_1d(FFT_REAL_SIZE, &instance->state->buffer21[0], instance->state->fft_buffer2, FFTW_FORWARD, FFTW_ESTIMATE);
    instance->state->fft_plan22 = fftw_plan_dft_1d(FFT_REAL_SIZE, &instance->state->buffer22[0], instance->state->fft_buffer2, FFTW_FORWARD, FFTW_ESTIMATE);
    
    instance->state->ifft_plan = fftw_plan_dft_1d(FFT_REAL_SIZE, instance->state->fft_buffer3, instance->state->fft_buffer4, FFTW_BACKWARD, FFTW_ESTIMATE);

    instance->state->convolution_buffer.resize(FFT_REAL_SIZE, 0);
    instance->state->convolution_buffer_head = FFT_REAL_SIZE - 1;

    return instance;
}

static void cleanup(plugin_t *instance) {
    fftw_destroy_plan (instance->state->fft_plan11);
    fftw_destroy_plan (instance->state->fft_plan12);
    fftw_destroy_plan (instance->state->fft_plan21);
    fftw_destroy_plan (instance->state->fft_plan22);
    fftw_destroy_plan (instance->state->ifft_plan);

    fftw_free (instance->state->spectrum1);
    fftw_free (instance->state->spectrum2);

    fftw_free (instance->state->buffer11);
    fftw_free (instance->state->buffer12);
    fftw_free (instance->state->buffer22);
    fftw_free (instance->state->buffer21);

    fftw_free (instance->state->response);

    fftw_free (instance->state->fft_buffer1);
    fftw_free (instance->state->fft_buffer2);
    fftw_free (instance->state->fft_buffer3);
    fftw_free (instance->state->fft_buffer4);
    // fftw_cleanup ();
    delete instance->state;
}

#define EPSILON 0.0001f

static void run(
    plugin_t *instance, uint32_t nframes, 
    const plugin_port_in_t in, const plugin_port_out_t out, 
    const plugin_port_analyze1_t analyze1, const plugin_port_analyze2_t analyze2, 
    const plugin_port_apply_t apply, const plugin_port_gain_t gain) 
{
    plugin_state_t *tinstance = instance->state;

    const bool previous_analyze1 = tinstance->previous_analyze1;
    const bool previous_analyze2 = tinstance->previous_analyze2;

    if (!previous_analyze1 && analyze1.data > 0) {
      // start collecting samples for spectrum 1
      std::cout << "start collecting spectrum1\n";
      memset(&tinstance->buffer11[0], 0, FFT_REAL_SIZE * 2);
      memset(&tinstance->buffer12[0], 0, FFT_REAL_SIZE * 2);
      memset(&tinstance->spectrum1[0], 0, FFT_COMPLEX_SIZE * 2);
      tinstance->buffer_head11 = 0;
      tinstance->buffer_head12 = FFT_REAL_SIZE/2;
      tinstance->number_of_ffts1 = 0;
    }

    if (!previous_analyze2 && analyze2.data > 0) {
      // start collecting samples for spectrum 2
      std::cout << "start collecting spectrum2\n";
      memset(&tinstance->buffer21[0], 0, FFT_REAL_SIZE * 2);
      memset(&tinstance->buffer22[0], 0, FFT_REAL_SIZE * 2);
      memset(&tinstance->spectrum2[0], 0, FFT_COMPLEX_SIZE * 2);
      tinstance->buffer_head21 = 0;
      tinstance->buffer_head22 = FFT_REAL_SIZE/2;
      tinstance->number_of_ffts2 = 0;
    }

    if ((previous_analyze1 && !(analyze1.data > 0)) || (previous_analyze2 && !(analyze2.data > 0))) {
      // calculate response
      std::cout << "calculate response\n";
      for (size_t index = 0; index < FFT_COMPLEX_SIZE; ++index) {
        tinstance->fft_buffer3[index][0] = (tinstance->spectrum2[index][0] / tinstance->spectrum1[index][0]) / FFT_REAL_SIZE;
        tinstance->fft_buffer3[index][1] = 0;

        tinstance->fft_buffer3[index][0] = log(tinstance->fft_buffer3[index][0]);
      }

      fftw_execute(tinstance->ifft_plan);

      // fold
      for (size_t index = 0; index < FFT_REAL_SIZE; ++index) {
        if (index == 0) {
          tinstance->buffer22[index][0] = tinstance->fft_buffer4[index][0];
          tinstance->buffer22[index][1] = 0;
        }
        if (index > 0 && index < FFT_REAL_SIZE/2) {
          tinstance->buffer22[index][0] = tinstance->fft_buffer4[index][0];
          tinstance->buffer22[index][1] = 0;
          tinstance->buffer22[index][0] += tinstance->fft_buffer4[FFT_REAL_SIZE - 1 - index][0];
          tinstance->buffer22[index][1] += 0;
        }
        if (index >= FFT_REAL_SIZE/2) {
          tinstance->buffer22[index][0] = 0;
          tinstance->buffer22[index][1] = 0;
        }
      }

      fftw_execute(tinstance->fft_plan22);

      for (size_t index = 0; index < FFT_REAL_SIZE; ++index) {
        tinstance->fft_buffer2[index][0] = exp(tinstance->fft_buffer2[index][0]);
        tinstance->fft_buffer2[index][1] = 0;
      }

      for (size_t index = 0; index < FFT_REAL_SIZE; ++index) {
          tinstance->fft_buffer3[index][0] = tinstance->fft_buffer2[index][0];
          tinstance->fft_buffer3[index][1] = tinstance->fft_buffer2[index][1];
      }

      fftw_execute (tinstance->ifft_plan);
      // ifftshift
      for (size_t index = 0; index < FFT_REAL_SIZE; ++index) tinstance->response[index][0] = tinstance->fft_buffer4[(index + FFT_REAL_SIZE/2) % FFT_REAL_SIZE][0];

      std::cout << "response: ";
      for (size_t index = 0; index < FFT_REAL_SIZE; ++index) std::cout << tinstance->response[index][0] << " ";
      std::cout << "\n";
      // for (size_t index = 0; index < FFT_REAL_SIZE/2; ++index) tinstance->response[index + FFT_REAL_SIZE/2][0] = 0;
    }

    for(uint32_t sample_index = 0; sample_index < nframes; ++sample_index) {
        if (analyze1.data > 0) {
            tinstance->buffer11[tinstance->buffer_head11][0] = in.data[sample_index] * tinstance->window[tinstance->buffer_head11];
            tinstance->buffer12[tinstance->buffer_head12][0] = in.data[sample_index] * tinstance->window[tinstance->buffer_head12];

            ++tinstance->buffer_head11;
            ++tinstance->buffer_head12;

            if (tinstance->buffer_head11 >= FFT_REAL_SIZE) {
                fftw_execute (tinstance->fft_plan11);
                for (size_t index = 0; index < FFT_COMPLEX_SIZE; ++index) {
                    tinstance->spectrum1[index][0] += sqrt(pow(tinstance->fft_buffer1[index][0], 2) + pow(tinstance->fft_buffer1[index][1], 2));
                    tinstance->spectrum1[index][1] = 0;
                }
                ++tinstance->number_of_ffts1;
                tinstance->buffer_head11 = 0;
            }

            if (tinstance->buffer_head12 >= FFT_REAL_SIZE) {
                fftw_execute (tinstance->fft_plan12);
                for (size_t index = 0; index < FFT_COMPLEX_SIZE; ++index) {
                    tinstance->spectrum1[index][0] += sqrt(pow(tinstance->fft_buffer1[index][0], 2) + pow(tinstance->fft_buffer1[index][1], 2));
                    tinstance->spectrum1[index][1] = 0;
                }
                ++tinstance->number_of_ffts1;
                tinstance->buffer_head12 = 0;
            }
        }

        if (analyze2.data > 0) {
            tinstance->buffer21[tinstance->buffer_head21][0] = in.data[sample_index] * tinstance->window[tinstance->buffer_head21];
            tinstance->buffer22[tinstance->buffer_head22][0] = in.data[sample_index] * tinstance->window[tinstance->buffer_head22];

            ++tinstance->buffer_head21;
            ++tinstance->buffer_head22;

            if (tinstance->buffer_head21 >= FFT_REAL_SIZE) {
                fftw_execute (tinstance->fft_plan21);
                for (size_t index = 0; index < FFT_COMPLEX_SIZE; ++index) {
                    tinstance->spectrum2[index][0] += sqrt(pow(tinstance->fft_buffer2[index][0], 2) + pow(tinstance->fft_buffer2[index][1], 2));
                    tinstance->spectrum2[index][1] = 0;
                }
                ++tinstance->number_of_ffts2;
                tinstance->buffer_head21 = 0;
            }

            if (tinstance->buffer_head22 >= FFT_REAL_SIZE) {
                fftw_execute (tinstance->fft_plan22);
                for (size_t index = 0; index < FFT_COMPLEX_SIZE; ++index) {
                    tinstance->spectrum2[index][0] += sqrt(pow(tinstance->fft_buffer2[index][0], 2) + pow(tinstance->fft_buffer2[index][1], 2));
                    tinstance->spectrum2[index][1] = 0;
                }
                ++tinstance->number_of_ffts2;
                tinstance->buffer_head22 = 0;
            }
        }

        tinstance->convolution_buffer[tinstance->convolution_buffer_head] = in.data[sample_index];
 
        out.data[sample_index] = in.data[sample_index];

        if (apply.data > 0) {
            out.data[sample_index] = 0;
            for (size_t index = 0; index < FFT_REAL_SIZE; ++index) {
                out.data[sample_index] += tinstance->response[index][0] * tinstance->convolution_buffer[(tinstance->convolution_buffer_head + index) % FFT_REAL_SIZE];
            }
        }        

        if (tinstance->convolution_buffer_head == 0) {
          tinstance->convolution_buffer_head = FFT_REAL_SIZE - 1;
        } else {
          --tinstance->convolution_buffer_head;
        }
    }

    tinstance->previous_analyze1 = analyze1.data > 0;
    tinstance->previous_analyze2 = analyze2.data > 0;
}

static const plugin_callbacks_t plugin_callbacks = {
    .instantiate = instantiate,
    .run = run,
    .cleanup = cleanup,
};

#include "generated/ttl2c_eq_match.c"

