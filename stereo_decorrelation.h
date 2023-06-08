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
        m_length (length)
    {
        init (1, 0);
    }

    void init (float decay, int seed)
    {
        m_left_response.resize (m_length);
        m_right_response.resize (m_length);
        
        m_left_response[0] = 0;
        m_right_response[0] = 0;
        
        std::mt19937 gen(seed); 
        std::normal_distribution<float> d(0, 0.1f);
        
        for (size_t index = 1; index < m_length; ++index)
        {
            float sample = d (gen);
            m_left_response[index] = sample * expf(-(index / decay));
            m_right_response[index] = -sample * expf(-(index / decay));
        }
    }
};

#endif