#ifndef EQ_MATCH_PLUGINS_EQ_MATCH_HH
#define EQ_MATCH_PLUGINS_EQ_MATCH_HH

#include <string.h>
#include <fftw3.h>
#include <cmath>
#include <complex>
#include <vector>
#include <iostream>

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

struct dft
{
    const size_t m_size;

    fftwf_complex *m_in;
    fftwf_complex *m_out;
    fftwf_plan m_plan_forward;
    fftwf_plan m_plan_backward;

    dft (const size_t size) :
        m_size (size)
    {
        m_in = fftwf_alloc_complex (m_size);
        m_out = fftwf_alloc_complex (m_size);

        m_plan_forward = fftwf_plan_dft_1d(m_size, m_in, m_out, FFTW_FORWARD, FFTW_ESTIMATE);
        m_plan_backward = fftwf_plan_dft_1d(m_size, m_in, m_out, FFTW_FORWARD, FFTW_ESTIMATE);
    }

    ~dft ()
    {
        fftwf_free (m_in);
        fftwf_free (m_out);
        fftwf_destroy_plan (m_plan_forward);
        fftwf_destroy_plan (m_plan_backward);
    }

    void fft (fftwf_complex *in, fftwf_complex *out, bool forward = true)
    {
        for (size_t index = 0; index < m_size; ++index)
        {
            m_in[index][0] = in[index][0];
            m_in[index][1] = in[index][1];
        }

        fftwf_execute (forward ? m_plan_forward : m_plan_backward);

        for (size_t index = 0; index < m_size; ++index)
        {
            out[index][0] = m_out[index][0];
            out[index][1] = m_out[index][1];
        }
    }

    void ifft (fftwf_complex *in, fftwf_complex *out)
    {
        fft (in, out, false);
    }

    // Convenience function
    void fft_from_real (EQ_MATCH_FLOAT *in, fftwf_complex *out)
    {
        for (size_t index = 0; index < m_size; ++index)
        {
            m_in[index][0] = in[index];
            m_in[index][1] = 0;
        }

        fftwf_execute (m_plan_forward);

        for (size_t index = 0; index < m_size; ++index)
        {
            out[index][0] = m_out[index][0];
            out[index][1] = m_out[index][1];
        }
    }
};

struct eq_match
{
    size_t m_fft_size;

    dft m_dft1;
    dft m_dft2;

    EQ_MATCH_FLOAT m_sample_rate;

    std::vector<EQ_MATCH_FLOAT> m_window;

    // Buffers to accumulate samples
    fftwf_complex *m_buffer11;
    int m_buffer_head11;
    fftwf_complex *m_buffer12;
    int m_buffer_head12;

    fftwf_complex *m_buffer21;
    int m_buffer_head21;
    fftwf_complex *m_buffer22;
    int m_buffer_head22;

    fftwf_complex *m_spectrum1;
    fftwf_complex *m_spectrum2;

    fftwf_plan m_fft_plan11;
    fftwf_plan m_fft_plan12;
    fftwf_plan m_fft_plan21;
    fftwf_plan m_fft_plan22;

    fftwf_plan m_fft_plan;
    fftwf_plan m_ifft_plan;

    fftwf_complex *m_fft_buffer1;
    size_t m_number_of_ffts1;

    fftwf_complex *m_fft_buffer2;
    size_t m_number_of_ffts2;

    fftwf_complex *m_fft_buffer3;
    fftwf_complex *m_fft_buffer4;
    fftwf_complex *m_fft_buffer5;

    EQ_MATCH_FLOAT *m_linear_phase_response;
    EQ_MATCH_FLOAT *m_minimum_phase_response;

    eq_match (size_t fft_size, EQ_MATCH_FLOAT sample_rate) :
        m_fft_size (fft_size),
        m_dft1 (fft_size),
        m_dft2 (4 * fft_size),
        m_sample_rate (sample_rate),
        m_window(fft_size, 0.f)
    {
        m_buffer11 = (fftwf_complex*)fftwf_alloc_complex (fft_size);
        m_buffer12 = (fftwf_complex*)fftwf_alloc_complex (fft_size);

        m_buffer21 = (fftwf_complex*)fftwf_alloc_complex (fft_size);
        m_buffer22 = (fftwf_complex*)fftwf_alloc_complex (fft_size);

        m_buffer_head11 = 0;
        m_buffer_head12 = fft_size/2;
        m_buffer_head21 = 0;
        m_buffer_head22 = fft_size/2;

        m_spectrum1 = (fftwf_complex*)fftwf_alloc_complex (fft_size);
        m_spectrum2 = (fftwf_complex*)fftwf_alloc_complex (fft_size);

        m_fft_buffer1 = (fftwf_complex*)fftwf_alloc_complex (fft_size);
        m_number_of_ffts1 = 0;

        m_fft_buffer2 = (fftwf_complex*)fftwf_alloc_complex (fft_size);
        m_number_of_ffts2 = 0;

        m_fft_buffer3 = (fftwf_complex*)fftwf_alloc_complex (fft_size);
        m_fft_buffer4 = (fftwf_complex*)fftwf_alloc_complex (fft_size);
        m_fft_buffer5 = (fftwf_complex*)fftwf_alloc_complex (fft_size);

        // hann window generation
        for (size_t index = 0; index < fft_size; ++index) {
            m_window[index] = pow(sin(M_PI * index / fft_size), 2);
        }

        m_fft_plan11 = fftwf_plan_dft_1d(fft_size, &m_buffer11[0], m_fft_buffer1, FFTW_FORWARD, FFTW_ESTIMATE);
        m_fft_plan12 = fftwf_plan_dft_1d(fft_size, &m_buffer12[0], m_fft_buffer1, FFTW_FORWARD, FFTW_ESTIMATE);

        m_fft_plan21 = fftwf_plan_dft_1d(fft_size, &m_buffer21[0], m_fft_buffer2, FFTW_FORWARD, FFTW_ESTIMATE);
        m_fft_plan22 = fftwf_plan_dft_1d(fft_size, &m_buffer22[0], m_fft_buffer2, FFTW_FORWARD, FFTW_ESTIMATE);

        m_fft_plan = fftwf_plan_dft_1d(fft_size, m_fft_buffer3, m_fft_buffer5, FFTW_FORWARD, FFTW_ESTIMATE);

        m_ifft_plan = fftwf_plan_dft_1d(fft_size, m_fft_buffer3, m_fft_buffer4, FFTW_BACKWARD, FFTW_ESTIMATE);

        m_linear_phase_response = (EQ_MATCH_FLOAT*)malloc (sizeof(EQ_MATCH_FLOAT) * fft_size);
        
        m_minimum_phase_response = (EQ_MATCH_FLOAT*)malloc (sizeof(EQ_MATCH_FLOAT) * fft_size);
        
        reset ();
    }

    ~eq_match ()
    {
        fftwf_destroy_plan (m_fft_plan11);
        fftwf_destroy_plan (m_fft_plan12);
        fftwf_destroy_plan (m_fft_plan21);
        fftwf_destroy_plan (m_fft_plan22);

        fftwf_destroy_plan (m_fft_plan);
        fftwf_destroy_plan (m_ifft_plan);

        fftwf_free (m_spectrum1);
        fftwf_free (m_spectrum2);

        fftwf_free (m_buffer11);
        fftwf_free (m_buffer12);
        fftwf_free (m_buffer22);
        fftwf_free (m_buffer21);

        fftwf_free (m_fft_buffer1);
        fftwf_free (m_fft_buffer2);
        fftwf_free (m_fft_buffer3);
        fftwf_free (m_fft_buffer4);
        fftwf_free (m_fft_buffer5);

        free (m_linear_phase_response);
        free (m_minimum_phase_response);
    }

    // Buffers up to FFT_SIZE samples and then calculates the spectrum and adds
    // it to spectrum1
    void add_frames_to_buffer1 (const EQ_MATCH_FLOAT *buffer, size_t number_of_samples)
    {
        add_frames_to_buffer (0, buffer, number_of_samples);
    }

    // Buffers up to FFT_SIZE samples and then calculates the spectrum and adds
    // it to spectrum2
    void add_frames_to_buffer2 (const EQ_MATCH_FLOAT *buffer, size_t number_of_samples)
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

        std::fill (m_linear_phase_response, m_linear_phase_response + m_fft_size, 0);
        std::fill (m_minimum_phase_response, m_minimum_phase_response + m_fft_size, 0);
    }

    void calculate_response ()
    {
        DBG_COMPLEX_VECTOR("spectrum1:", m_spectrum1, m_fft_size)
        DBG_COMPLEX_VECTOR("spectrum2:", m_spectrum2, m_fft_size)

        for (size_t index = 0; index < m_fft_size; ++index) {
            m_fft_buffer3[index][0] = ((m_spectrum2[index][0] / m_number_of_ffts2) / (m_spectrum1[index][0] / m_number_of_ffts1));
            m_fft_buffer3[index][1] = 0;
        }

        DBG_COMPLEX_VECTOR("ratio:", m_fft_buffer3, m_fft_size)

        fftwf_execute(m_ifft_plan);
        for (size_t index = 0; index < m_fft_size; ++index)
        {
            m_fft_buffer4[index][0] /= m_fft_size;
        }

        // IFFT-shift
        for (size_t index = 0; index < m_fft_size; ++index)
        {
            m_linear_phase_response[(index + (m_fft_size/2)) % m_fft_size] = m_fft_buffer4[index][0];
        }

        // exp(fft(fold(ifft(log(s)))))
        for (size_t index = 0; index < m_fft_size; ++index) {
            m_fft_buffer3[index][0] = log(m_fft_buffer3[index][0]);
        }

        DBG_COMPLEX_VECTOR("log:", m_fft_buffer3, m_fft_size)

        // 3 -> 4
        fftwf_execute(m_ifft_plan);

        for (size_t index = 0; index < m_fft_size; ++index)
        {
            m_fft_buffer4[index][0] /= m_fft_size;
        }

        DBG_COMPLEX_VECTOR("ifft:", m_fft_buffer4, m_fft_size)

        // fold
        for (size_t index = 0; index < m_fft_size; ++index) {
            if (index == 0) {
                m_fft_buffer3[index][0] = m_fft_buffer4[index][0];
                m_fft_buffer3[index][1] = 0;
            }
            if (index > 0 && index < m_fft_size/2) {
                m_fft_buffer3[index][0] = m_fft_buffer4[index][0];
                m_fft_buffer3[index][1] = 0;
                m_fft_buffer3[index][0] += m_fft_buffer4[m_fft_size - index][0];
            }
            if (index == m_fft_size/2) {
                m_fft_buffer3[index][0] = m_fft_buffer4[index][0];
            }
            if (index > m_fft_size/2) {
                m_fft_buffer3[index][0] = 0;
                m_fft_buffer3[index][1] = 0;
            }
        }

        DBG_COMPLEX_VECTOR("fold:", m_fft_buffer3, m_fft_size)

        // 3 -> 5
        fftwf_execute (m_fft_plan);

        for (size_t index = 0; index < m_fft_size; ++index) {
            std::complex<EQ_MATCH_FLOAT> c(m_fft_buffer5[index][0], m_fft_buffer5[index][1]);
            std::complex<EQ_MATCH_FLOAT> c2 = std::exp(c);
            m_fft_buffer3[index][0] = std::real(c2);
            m_fft_buffer3[index][1] = std::imag(c2);
        }

        DBG_COMPLEX_VECTOR("exp:", m_fft_buffer3, m_fft_size)

        // 3 -> 4
        fftwf_execute (m_ifft_plan);
        for (size_t index = 0; index < m_fft_size; ++index) m_fft_buffer4[index][0] /= m_fft_size;

        for (size_t index = 0; index < m_fft_size; ++index) m_minimum_phase_response[index] = m_fft_buffer4[index][0];

        DBG_REAL_VECTOR("linear phase response:", m_linear_phase_response, m_fft_size)
        
        DBG_REAL_VECTOR("minimum phase response:", m_minimum_phase_response, m_fft_size)
    }

protected:
    void add_frames_to_buffer (size_t buffer_index, const EQ_MATCH_FLOAT *buffer, size_t number_of_samples)
    {
        fftwf_complex *buffer1 = (0 == buffer_index) ? m_buffer11 : m_buffer21;
        fftwf_complex *buffer2 = (0 == buffer_index) ? m_buffer12 : m_buffer22;

        fftwf_plan fft_plan1 = (0 == buffer_index) ? m_fft_plan11 : m_fft_plan21;
        fftwf_plan fft_plan2 = (0 == buffer_index) ? m_fft_plan12 : m_fft_plan22;

        fftwf_complex *fft_buffer = (0 == buffer_index) ? m_fft_buffer1 : m_fft_buffer2;

        fftwf_complex *spectrum = (0 == buffer_index) ? m_spectrum1 : m_spectrum2;

        size_t &number_of_ffts = (0 == buffer_index) ? m_number_of_ffts1 : m_number_of_ffts2;

        int &buffer_head1 = 0 == (buffer_index) ? m_buffer_head11 : m_buffer_head21;
        int &buffer_head2 = 0 == (buffer_index) ? m_buffer_head12 : m_buffer_head22;

        for (size_t sample_index = 0; sample_index < number_of_samples; ++sample_index)
        {
            buffer1[buffer_head1][0] = buffer[sample_index] * m_window[buffer_head1];
            buffer2[buffer_head2][0] = buffer[sample_index] * m_window[buffer_head2];

            ++buffer_head1;
            ++buffer_head2;

            if (buffer_head1 >= (int)m_fft_size) {
                    // DBG("execute " << buffer_index << "-1\n")
                    fftwf_execute (fft_plan1);
                    for (size_t index = 0; index < m_fft_size; ++index) {
                            spectrum[index][0] += sqrt(pow(fft_buffer[index][0], 2) + pow(fft_buffer[index][1], 2));
                            spectrum[index][1] = 0;
                    }
                    ++number_of_ffts;
                    buffer_head1 = 0;
            }

            if (buffer_head2 >= (int)m_fft_size) {
                    // DBG("execute " << buffer_index << "-2\n")
                    fftwf_execute (fft_plan2);
                    for (size_t index = 0; index < m_fft_size; ++index) {
                            spectrum[index][0] += sqrt(pow(fft_buffer[index][0], 2) + pow(fft_buffer[index][1], 2));
                            spectrum[index][1] = 0;
                    }
                    ++number_of_ffts;
                    buffer_head2 = 0;
            }
        }
    }

    void reset_buffer (size_t buffer_index)
    {
        fftwf_complex *buffer1 = (0 == buffer_index) ? m_buffer11 : m_buffer21;
        fftwf_complex *buffer2 = (0 == buffer_index) ? m_buffer12 : m_buffer22;

        fftwf_complex *spectrum = (0 == buffer_index) ? m_spectrum1 : m_spectrum2;

        size_t &number_of_ffts = (0 == buffer_index) ? m_number_of_ffts1 : m_number_of_ffts2;

        int &buffer_head1 = 0 == (buffer_index) ? m_buffer_head11 : m_buffer_head21;
        int &buffer_head2 = 0 == (buffer_index) ? m_buffer_head12 : m_buffer_head22;

        memset(buffer1, 0, sizeof(EQ_MATCH_FLOAT) * m_fft_size * 2);
        memset(buffer2, 0, sizeof(EQ_MATCH_FLOAT) * m_fft_size * 2);
        memset(spectrum, 0, sizeof(EQ_MATCH_FLOAT) * m_fft_size * 2);

        buffer_head1 = 0;
        buffer_head2 = m_fft_size/2;

        number_of_ffts = 0;
    }
};

#undef fftwf_complex
#undef fftwf_plan
#undef fftwf_plan_dft_1d
#undef fftwf_alloc_complex
#undef fftwf_free
#undef fftwf_destroy_plan
#undef fftwf_execute

#endif
