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

    for (size_t index = 0; index < length; ++index)
    {
        float sample = d (gen);
        response_out[index] = sample * exp(-(index / decay));
    }

    if (whiten)
    {
        whiten_by_minimum_phase_spectrum (response_out, response_out);
    }

    float dc = 0;
    for (size_t index = 0; index < length; ++index)
    {
        dc += response_out[index];
    }
    dc /= length;

    float power = 0;
    for (size_t index = 0; index < length; ++index)
    {
        // response_out[index] -= dc;
        power += response_out[index] * response_out[index];
    }
    power = sqrtf (power);

    for (size_t index = 0; index < length; ++index)
    {
        response_out[index] /= power;
    }
}

void create_exponential_white_noise_burst_stereo_pair
(
    const int random_seed_left,
    const int random_seed_right,
    const float decay,
    const bool whiten,
    const bool sum_to_mono,
    std::vector<float> &response_left_out,
    std::vector<float> &response_right_out
)
{
    assert (response_left_out.size () == response_right_out.size ());
    const size_t length = response_left_out.size ();
    assert (length >= 8);

    create_exponential_white_noise_burst (random_seed_left, decay, whiten, response_left_out);

    if (sum_to_mono)
    {
        for (size_t index = 0; index < length; ++index)
        {
            response_right_out[index] = -response_left_out[index];
        }
    }
    else
    {
        create_exponential_white_noise_burst (random_seed_right, decay, whiten, response_right_out);
    }
}


#endif
