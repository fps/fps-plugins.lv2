#ifndef EQ_MATCH_PLUGINS_DECORRELATION_HH
#define EQ_MATCH_PLUGINS_DECORRELATION_HH

#include <string.h>
#include <vector>
#include <random>
#include <cmath>
#include <iostream>

#include "minimum_phase_spectrum.h"
#include "dft.h"

void whiten_by_minimum_phase_spectrum
(
    const std::vector<float> &in,
    std::vector<float> &out
)
{
    assert (in.size () == out.size ());
    assert (in.size () >= 8);

    size_t fft_size = in.size ();

    dft the_dft (fft_size);
    dft_buffer buffer1 (fft_size);
    dft_buffer buffer2 (fft_size);
    dft_buffer buffer3 (fft_size);

    for (size_t index = 0; index < fft_size; ++index)
    {
        buffer1.m[index][0] = in[index];
        buffer1.m[index][1] = 0;
    }

    the_dft.fft (buffer1, buffer2);

    for (size_t index = 0; index < fft_size; ++index)
    {
        buffer2.m[index][0] = sqrt(pow(buffer2.m[index][0], 2) + pow(buffer2.m[index][1], 2));
        buffer2.m[index][1] = 0;
    }

    create_minimum_phase_spectrum (buffer2, buffer3);

    for (size_t index = 0; index < fft_size; ++index)
    {
        buffer1.m[index][0] = in[index];
        buffer1.m[index][1] = 0;
    }

    the_dft.fft (buffer1, buffer2);

    for (size_t index = 0; index < fft_size; ++index)
    {
        std::complex c1 (buffer2.m[index][0], buffer2.m[index][1]);
        std::complex c2 (buffer3.m[index][0], buffer3.m[index][1]);

        std::complex ratio = c1 / c2;

        buffer1.m[index][0] = std::real (ratio);
        buffer1.m[index][1] = std::imag (ratio);
    }

    the_dft.ifft (buffer1, buffer2);

    for (size_t index = 0; index < fft_size; ++index)
    {
        out[index] = buffer2.m[index][0];
    }
}

void create_exponential_white_noise_burst
(
    const int random_seed,
    const float decay,
    const bool whiten,
    std::vector<float> &response_out
)
{
    const size_t length = response_out.size ();
    assert (length >= 8);

    std::mt19937 gen(random_seed);
    std::normal_distribution<float> d(0, 1.0f);

    float dc = 0;
    for (size_t index = 0; index < length; ++index)
    {
        float sample = d (gen);
        response_out[index] = sample * exp(-(index / decay));
        dc += sample;
    }

    float power = 0;
    for (size_t index = 0; index < length; ++index)
    {
        // response_out[index] -= dc / length;
        power += response_out[index] * response_out[index];
    }
    power = sqrtf (power);

    for (size_t index = 0; index < length; ++index)
    {
        // response_out[index] /= power;
    }

    if (whiten)
    {
        whiten_by_minimum_phase_spectrum (response_out, response_out);
    }
}

struct stereo_decorrelation
{
    size_t m_length;

    std::vector<float> m_left_response;
    std::vector<float> m_right_response;
    
    stereo_decorrelation (size_t length) :
        m_length (length),
        m_left_response (length, 0),
        m_right_response (length, 0)
    {
        init (1, 11, 12, true, false);
    }

    void init (float decay, int seed1, int seed2, bool whitening, bool sum_to_mono)
    {
        m_left_response.resize (m_length);
        m_right_response.resize (m_length);
        
        std::mt19937 gen1(seed1);
        std::mt19937 gen2(seed2);

        std::normal_distribution<float> d(0, 0.5f);
        
        float dc1 = 0;
        float dc2 = 0;
        for (size_t index = 0; index < m_length; ++index)
        {
            float sample1 = d (gen1);
            float sample2 = d (gen2);
            m_left_response[index] = sample1 * expf(-(index / decay));
            dc1 += m_left_response[index];
            m_right_response[index] = sample2 * expf(-(index / decay));
            dc2 += m_right_response[index];
        }

        // remove DC
        float power1 = 0;
        float power2 = 0;
        for (size_t index = 0; index < m_length; ++index)
        {
            m_left_response[index] -= dc1 / m_length;
            power1 += powf(m_left_response[index], 2);
            m_right_response[index] -= dc1 / m_length;
            power1 += powf(m_right_response[index], 2);
        }
        power1 = sqrtf(power1);
        power2 = sqrtf(power2);

        // normalize to unit power
        for (size_t index = 0; index < m_length; ++index)
        {
            m_left_response[index] /= power1;
            m_right_response[index] /= power2;
        }

        if (sum_to_mono)
        {
            // invert right channel
            for (size_t index = 0; index < m_length; ++index)
            {
                m_right_response[index] = -m_left_response[index];
            }
        }
    }
};

#endif
