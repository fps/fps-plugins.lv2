#ifndef EQ_MATCH_PLUGINS_EQ_MATCH_HH
#define EQ_MATCH_PLUGINS_EQ_MATCH_HH

#include <string.h>
#include <vector>
#include <random>
#include <cmath>
#include <iostream>

#include "dft.h"

/**
 * This whole class is not realtime safe!
 */


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
        init (1, 0);
    }

    void init (float decay, int seed)
    {
        m_left_response.resize (m_length);
        m_right_response.resize (m_length);
        
        std::mt19937 gen(seed);
        std::normal_distribution<float> d(0, 0.5f);
        
        float dc = 0;
        for (size_t index = 0; index < m_length; ++index)
        {
            float sample = d (gen);
            m_left_response[index] = sample * expf(-(index / decay));
            dc += m_left_response[index];
        }

        // remove DC
        float power = 0;
        for (size_t index = 0; index < m_length; ++index)
        {
            m_left_response[index] -= dc / m_length;
            power += powf(m_left_response[index], 2);
        }
        power = sqrtf(power);

        // normalize to unit power
        for (size_t index = 0; index < m_length; ++index)
        {
            m_left_response[index] /= power;
        }

        // invert right channel
        for (size_t index = 0; index < m_length; ++index)
        {
            m_right_response[index] = -m_left_response[index];
        }
    }
};

#endif
