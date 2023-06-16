#ifndef EQ_MATCH_PLUGINS_DECORRELATION_HH
#define EQ_MATCH_PLUGINS_DECORRELATION_HH

#include <string.h>
#include <vector>
#include <random>
#include <cmath>
#include <iostream>

#include "dft.h"

/**
 * This whole class is not realtime safe!
 */

struct exponential_white_noise_decorrelator
{
    size_t m_length;
    std::vector<float> m_response;

    exponential_white_noise_decorrelator (size_t length) :
        m_length (length),
        m_response (length, 0)
    {

    }

    void init (bool whiten, float decay, int random_seed)
    {
        std::mt19937 gen(random_seed);
        std::normal_distribution<float> d(0, 1.0f);

        float dc = 0;
        for (size_t index = 0; index < m_length; ++index)
        {
            float sample = d (gen);
            m_response[index] = sample * expf(-index / decay);
            dc += sample;
        }

        float power = 0;
        for (size_t index = 0; index < m_length; ++index)
        {
            m_response[index] -= dc / m_length;
            power += m_response[index] * m_response[index];
        }
        power = sqrtf (power);

        for (size_t index = 0; index < m_length; ++index)
        {
            m_response[index] /= power;
        }

        if (whiten)
        {

        }
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
