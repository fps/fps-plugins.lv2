#ifndef EQ_MATCH_PLUGINS_EQ_MATCH_HH
#define EQ_MATCH_PLUGINS_EQ_MATCH_HH

#include <string.h>
#include <fftw3.h>
#include <cmath>
#include <complex>
#include <vector>
#include <iostream>

#include "dft.h"
#include "minimum_phase_spectrum.h"

/**
 * This whole class is not realtime safe!
 *
 * How to:
 * (1) Construct an object of type eq_match
 * (2) Call reset_buffer1 () to reset the spectrum estimation buffer1
 * (3) Call add_frames_to_buffer1 consecutively to estimate spectrum1
 * (4) Call reset_buffer2 () to reset the spectrum estimation buffer2
 * (5) Call add_frames_to_buffer2 consecutively to estimate spectrum2
 * (6) Call calculate_response ()
 * (7) Get the calculated IRs from m_response (linear phase) or
 *         m_minimum_phase_response (which is minimum phase)
 */

#define EQ_MATCH_MINIMUM_AMPLITUDE 0.0000001

struct eq_match
{
    size_t m_fft_size;
    size_t m_extended_fft_size;

    dft m_dft1;
    dft m_dft2;

    std::vector<float> m_window;

    // Buffers to accumulate samples
    dft_buffer m_spectrum_buffer11;
    int m_spectrum_buffer_head11;
    dft_buffer m_spectrum_buffer12;
    int m_spectrum_buffer_head12;

    dft_buffer m_spectrum_buffer21;
    int m_spectrum_buffer_head21;
    dft_buffer m_spectrum_buffer22;
    int m_spectrum_buffer_head22;

    dft_buffer m_average_spectrum1;
    dft_buffer m_average_spectrum2;

    std::vector<float> m_linear_phase_response;
    std::vector<float> m_minimum_phase_response;

    eq_match (size_t fft_size) :
        m_fft_size (fft_size),
        m_extended_fft_size (fft_size * 4),
        m_dft1 (fft_size),
        m_dft2 (4 * fft_size),
        m_window(fft_size, 0.f),

        m_spectrum_buffer11 (fft_size),
        m_spectrum_buffer12 (fft_size),
        m_spectrum_buffer21 (fft_size),
        m_spectrum_buffer22 (fft_size),

        m_average_spectrum1 (fft_size),
        m_average_spectrum2 (fft_size),

        m_linear_phase_response (fft_size),
        m_minimum_phase_response (fft_size)
    {
        // window generation
        // const float a0 = 25.0f/46.0f;
        for (size_t index = 0; index < fft_size; ++index) {
            m_window[index] = pow(sin(M_PI * index / fft_size), 2);
            // m_window[index] = a0 - (1.0f - a0) * cos(2*M_PI*index/fft_size);
        }

        DBG_REAL_VECTOR("window:", m_window, m_fft_size)
        reset ();
    }

    ~eq_match ()
    {

    }

    // Buffers up to FFT_SIZE samples and then calculates the spectrum and adds
    // it to spectrum1
    void add_frames_to_buffer1 (const float *buffer, size_t number_of_samples)
    {
        add_frames_to_buffer (0, buffer, number_of_samples);
    }

    // Buffers up to FFT_SIZE samples and then calculates the spectrum and adds
    // it to spectrum2
    void add_frames_to_buffer2 (const float *buffer, size_t number_of_samples)
    {
        add_frames_to_buffer (1, buffer, number_of_samples);
    }

    void reset_buffer1 ()
    {
        reset_buffer (0);
    }

    void reset_buffer2 ()
    {
        reset_buffer (1);
    }

    void reset ()
    {
        reset_buffer1 ();
        reset_buffer2 ();

        std::fill (m_linear_phase_response.begin (), m_linear_phase_response.end (), 0);
        std::fill (m_minimum_phase_response.begin (), m_minimum_phase_response.end (), 0);
    }

    void calculate_response ()
    {
        DBG_COMPLEX_VECTOR("spectrum1:", m_average_spectrum1.m, m_fft_size)
        DBG_COMPLEX_VECTOR("spectrum2:", m_average_spectrum2.m, m_fft_size)

        float max_amplitude1 = 0;
        for (size_t index = 0; index < m_fft_size; ++index)
        {
            if (m_average_spectrum1.m[index][0] < EQ_MATCH_MINIMUM_AMPLITUDE)
            {
                m_average_spectrum1.m[index][0] = EQ_MATCH_MINIMUM_AMPLITUDE;
            }
            // power1 += m_average_spectrum1.m[index][0] * m_average_spectrum1.m[index][0];
            if (m_average_spectrum1.m[index][0] > max_amplitude1) max_amplitude1 = m_average_spectrum1.m[index][0];
        }
        
        float max_amplitude2 = 0;
        for (size_t index = 0; index < m_fft_size; ++index)
        {
            if (m_average_spectrum2.m[index][0] < EQ_MATCH_MINIMUM_AMPLITUDE)
            {
                m_average_spectrum2.m[index][0] = EQ_MATCH_MINIMUM_AMPLITUDE;
            }
            // power2 += m_average_spectrum2.m[index][0] * m_average_spectrum2.m[index][0];
            if (m_average_spectrum2.m[index][0] > max_amplitude2) max_amplitude2 = m_average_spectrum2.m[index][0];
        }

        dft_buffer spectrum_ratio (m_fft_size);
        for (size_t index = 0; index < m_fft_size; ++index) {
            spectrum_ratio.m[index][0] =
              (m_average_spectrum2.m[index][0] / max_amplitude1) / (m_average_spectrum1.m[index][0] / max_amplitude2);

            spectrum_ratio.m[index][1] = 0;
        }

        DBG_COMPLEX_VECTOR("ratio:", spectrum_ratio.m, m_fft_size)

        dft_buffer ifft_spectrum_ratio (m_fft_size);
        m_dft1.ifft (spectrum_ratio, ifft_spectrum_ratio);

        // IFFT-shift
        for (size_t index = 0; index < m_fft_size; ++index)
        {
            m_linear_phase_response[(index + (m_fft_size/2)) % m_fft_size] = ifft_spectrum_ratio.m[index][0];
        }

        create_minimum_phase_response (m_linear_phase_response, m_minimum_phase_response);

        DBG_REAL_VECTOR("linear phase response:", m_linear_phase_response, m_fft_size)
        
        DBG_REAL_VECTOR("minimum phase response:", m_minimum_phase_response, m_fft_size)
    }

protected:
    void add_frames_to_buffer (size_t buffer_index, const float *buffer, size_t number_of_samples)
    {
        dft_buffer &buffer1 = (0 == buffer_index) ? m_spectrum_buffer11 : m_spectrum_buffer21;
        dft_buffer &buffer2 = (0 == buffer_index) ? m_spectrum_buffer12 : m_spectrum_buffer22;

        dft_buffer &spectrum = (0 == buffer_index) ? m_average_spectrum1 : m_average_spectrum2;

        int &buffer_head1 = 0 == (buffer_index) ? m_spectrum_buffer_head11 : m_spectrum_buffer_head21;
        int &buffer_head2 = 0 == (buffer_index) ? m_spectrum_buffer_head12 : m_spectrum_buffer_head22;

        for (size_t sample_index = 0; sample_index < number_of_samples; ++sample_index)
        {
            buffer1.m[buffer_head1][0] = buffer[sample_index] * m_window[buffer_head1];
            buffer2.m[buffer_head2][0] = buffer[sample_index] * m_window[buffer_head2];

            ++buffer_head1;
            ++buffer_head2;

            if (buffer_head1 >= (int)m_fft_size) {
                    // DBG("execute " << buffer_index << "-1\n")
                    dft_buffer current_spectrum (m_fft_size);
                    m_dft1.fft (buffer1, current_spectrum);
                    for (size_t index = 0; index < m_fft_size; ++index) {
                            spectrum.m[index][0] += sqrtf(powf(current_spectrum.m[index][0], 2) + powf(current_spectrum.m[index][1], 2));
                            spectrum.m[index][1] = 0;
                    }
                    buffer_head1 = 0;
            }

            if (buffer_head2 >= (int)m_fft_size) {
                    // DBG("execute " << buffer_index << "-2\n")
                    dft_buffer current_spectrum (m_fft_size);
                    m_dft1.fft (buffer2, current_spectrum);
                    for (size_t index = 0; index < m_fft_size; ++index) {
                            spectrum.m[index][0] += sqrtf(powf(current_spectrum.m[index][0], 2) + powf(current_spectrum.m[index][1], 2));
                            spectrum.m[index][1] = 0;
                    }
                    buffer_head2 = 0;
            }
        }
    }

    void reset_buffer (size_t buffer_index)
    {
        dft_buffer &buffer1 = (0 == buffer_index) ? m_spectrum_buffer11 : m_spectrum_buffer21;
        dft_buffer &buffer2 = (0 == buffer_index) ? m_spectrum_buffer12 : m_spectrum_buffer22;

        dft_buffer &spectrum = (0 == buffer_index) ? m_average_spectrum1 : m_average_spectrum2;

        int &buffer_head1 = 0 == (buffer_index) ? m_spectrum_buffer_head11 : m_spectrum_buffer_head21;
        int &buffer_head2 = 0 == (buffer_index) ? m_spectrum_buffer_head12 : m_spectrum_buffer_head22;

        memset(buffer1.m, 0, sizeof(float) * m_fft_size * 2);
        memset(buffer2.m, 0, sizeof(float) * m_fft_size * 2);
        memset(spectrum.m, 0, sizeof(float) * m_fft_size * 2);

        buffer_head1 = 0;
        buffer_head2 = m_fft_size/2;
    }
};

#endif
