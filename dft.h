#ifndef FPS_PLUGINS_DF_BUFFER_HH
#define FPS_PLUGINS_DF_BUFFER_HH

#include <fftw3.h>
#include <iostream>
#include <cmath>
#include <complex>

#ifdef DEBUG
    #define DBG(x) \
    {\
        std::cerr << x;\
    }

    #define DBG_REAL_VECTOR(prefix, x, len) \
    { \
        std::cerr << "############ " << prefix << "\n"; \
        for (size_t index = 0; index < len; ++index) \
        { \
            std::cerr << x[index] << " "; \
        } \
        std::cerr << "\n";\
    }

    #define DBG_COMPLEX_VECTOR(prefix, x, len) \
        { \
            std::cerr << "############ " << prefix << " real:\n"; \
            for (size_t index = 0; index < len; ++index) \
            { \
                std::cerr << x[index][0] << " ";\
            } \
            std::cerr << "\n"; \
            std::cerr << "############ " << prefix << " imag:\n"; \
            for (size_t index = 0; index < len; ++index) \
            { \
                std::cerr << x[index][1] << " "; \
            } \
            std::cerr << "\n";\
        }
#else
    #define DBG(x) {}

    #define DBG_REAL_VECTOR(prefix, x, len) {}

    #define DBG_COMPLEX_VECTOR(prefix, x, len) {}
#endif

struct dft_buffer
{
    size_t m_size;
    fftwf_complex *m;

    dft_buffer (size_t size) :
        m_size (size),
        m(fftwf_alloc_complex (size))
    {

    }

    ~dft_buffer ()
    {
      fftwf_free (m);
    }

private:
    dft_buffer (const dft_buffer &);
    dft_buffer &operator= (const dft_buffer &);
};

struct dft
{
    const size_t m_size;

    dft_buffer m_in;
    dft_buffer m_out;
    fftwf_plan m_plan_forward;
    fftwf_plan m_plan_backward;

    dft (const size_t size) :
        m_size (size),
        m_in (size),
        m_out (size)
    {
        m_plan_forward = fftwf_plan_dft_1d(m_size, m_in.m, m_out.m, FFTW_FORWARD, FFTW_ESTIMATE);
        m_plan_backward = fftwf_plan_dft_1d(m_size, m_in.m, m_out.m, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    ~dft ()
    {
        fftwf_destroy_plan (m_plan_forward);
        fftwf_destroy_plan (m_plan_backward);
    }

    void fft (const dft_buffer &in, dft_buffer &out, bool forward = true)
    {
        assert (in.m_size == m_size);
        assert (out.m_size == m_size);

        for (size_t index = 0; index < m_size; ++index)
        {
            m_in.m[index][0] = in.m[index][0];
            m_in.m[index][1] = in.m[index][1];
        }

        fftwf_execute (forward ? m_plan_forward : m_plan_backward);

        for (size_t index = 0; index < m_size; ++index)
        {
            out.m[index][0] = m_out.m[index][0];
            out.m[index][1] = m_out.m[index][1];
        }
    }

    void ifft (const dft_buffer &in, dft_buffer &out)
    {
        assert (in.m_size == m_size);
        assert (out.m_size == m_size);

        fft (in, out, false);

        // normalize the ifft path
        for (size_t index = 0; index < m_size; ++index)
        {
            out.m[index][0] /= m_size;
            out.m[index][1] /= m_size;
        }
    }

    // Convenience function
    void real_fft (const std::vector<float> &in, dft_buffer &out)
    {
        assert (in.size () == out.m_size);

        for (size_t index = 0; index < m_size; ++index)
        {
            m_in.m[index][0] = in[index];
            m_in.m[index][1] = 0;
        }

        fftwf_execute (m_plan_forward);

        for (size_t index = 0; index < m_size; ++index)
        {
            out.m[index][0] = m_out.m[index][0];
            out.m[index][1] = m_out.m[index][1];
        }
    }
};


#endif
