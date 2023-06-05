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
        m_plan_backward = fftwf_plan_dft_1d(m_size, m_in.m, m_out.m, FFTW_FORWARD, FFTW_ESTIMATE);
    }

    ~dft ()
    {
        fftwf_destroy_plan (m_plan_forward);
        fftwf_destroy_plan (m_plan_backward);
    }

    void fft (const dft_buffer &in, dft_buffer &out, bool forward = true)
    {
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
        fft (in, out, false);
    }

    // Convenience function
    void fft_from_real (float *in, fftwf_complex *out)
    {
        for (size_t index = 0; index < m_size; ++index)
        {
            m_in.m[index][0] = in[index];
            m_in.m[index][1] = 0;
        }

        fftwf_execute (m_plan_forward);

        for (size_t index = 0; index < m_size; ++index)
        {
            out[index][0] = m_out.m[index][0];
            out[index][1] = m_out.m[index][1];
        }
    }
};

struct eq_match
{
    size_t m_fft_size;

    dft m_dft1;
    dft m_dft2;

    float m_sample_rate;

    std::vector<float> m_window;

    // Buffers to accumulate samples
    dft_buffer m_buffer11;
    int m_buffer_head11;
    dft_buffer m_buffer12;
    int m_buffer_head12;

    dft_buffer m_buffer21;
    int m_buffer_head21;
    dft_buffer m_buffer22;
    int m_buffer_head22;

    dft_buffer m_spectrum1;
    dft_buffer m_spectrum2;

    dft_buffer m_spectrum_ratio;
    dft_buffer m_ifft_spectrum_ratio;
    dft_buffer m_log;
    dft_buffer m_ifft_log;
    dft_buffer m_fold_ifft_log;
    dft_buffer m_fft_fold_ifft_log;
    dft_buffer m_exp_fft_fold_ifft_log;
    dft_buffer m_ifft_exp_fft_fold_ifft_log;

    dft_buffer m_fft_buffer1;
    size_t m_number_of_ffts1;

    dft_buffer m_fft_buffer2;
    size_t m_number_of_ffts2;

#if 0
    fftwf_plan m_fft_plan11;
    fftwf_plan m_fft_plan12;
    fftwf_plan m_fft_plan21;
    fftwf_plan m_fft_plan22;

    fftwf_plan m_fft_plan;
    fftwf_plan m_ifft_plan;


    fftwf_complex *m_fft_buffer3;
    fftwf_complex *m_fft_buffer4;
    fftwf_complex *m_fft_buffer5;
#endif

    std::vector<float> m_linear_phase_response;
    std::vector<float> m_minimum_phase_response;

    eq_match (size_t fft_size, float sample_rate) :
        m_fft_size (fft_size),
        m_dft1 (fft_size),
        m_dft2 (4 * fft_size),
        m_sample_rate (sample_rate),
        m_window(fft_size, 0.f),
        m_buffer11 (fft_size),
        m_buffer12 (fft_size),
        m_buffer21 (fft_size),
        m_buffer22 (fft_size),
        m_spectrum1 (fft_size),
        m_spectrum2 (fft_size),
        m_spectrum_ratio (fft_size),
        m_ifft_spectrum_ratio (fft_size),
        m_log (fft_size),
        m_ifft_log (fft_size),
        m_fold_ifft_log (fft_size),
        m_fft_fold_ifft_log (fft_size),
        m_exp_fft_fold_ifft_log (fft_size),
        m_ifft_exp_fft_fold_ifft_log (fft_size),
        m_fft_buffer1 (fft_size),
        m_fft_buffer2 (fft_size),
        m_linear_phase_response (fft_size),
        m_minimum_phase_response (fft_size)
    {
        // hann window generation
        for (size_t index = 0; index < fft_size; ++index) {
            m_window[index] = pow(sin(M_PI * index / fft_size), 2);
        }

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
        DBG_COMPLEX_VECTOR("spectrum1:", m_spectrum1, m_fft_size)
        DBG_COMPLEX_VECTOR("spectrum2:", m_spectrum2, m_fft_size)

        for (size_t index = 0; index < m_fft_size; ++index) {
            m_spectrum_ratio.m[index][0] = ((m_spectrum2.m[index][0] / m_number_of_ffts2) / (m_spectrum1.m[index][0] / m_number_of_ffts1));
            m_spectrum_ratio.m[index][1] = 0;
        }

        DBG_COMPLEX_VECTOR("ratio:", m_spectrum_ratio.m, m_fft_size)

        m_dft1.fft (m_spectrum_ratio, m_ifft_spectrum_ratio);

        for (size_t index = 0; index < m_fft_size; ++index)
        {
            m_ifft_spectrum_ratio.m[index][0] /= m_fft_size;
        }

        // IFFT-shift
        for (size_t index = 0; index < m_fft_size; ++index)
        {
            m_linear_phase_response[(index + (m_fft_size/2)) % m_fft_size] = m_ifft_spectrum_ratio.m[index][0];
        }

        // exp(fft(fold(ifft(log(s)))))

        for (size_t index = 0; index < m_fft_size; ++index) {
            m_log.m[index][0] = log(m_spectrum_ratio.m[index][0]);
        }

        DBG_COMPLEX_VECTOR("log:", m_spectrum_ratio.m, m_fft_size)

        // 3 -> 4
        m_dft1.ifft (m_log, m_ifft_log);

        for (size_t index = 0; index < m_fft_size; ++index)
        {
            m_ifft_log.m[index][0] /= m_fft_size;
        }

        DBG_COMPLEX_VECTOR("ifft:", m_ifft_log, m_fft_size)

        // fold
        for (size_t index = 0; index < m_fft_size; ++index) {
            if (index == 0) {
                m_fold_ifft_log.m[index][0] = m_ifft_log.m[index][0];
                m_fold_ifft_log.m[index][1] = 0;
            }
            if (index > 0 && index < m_fft_size/2) {
                m_fold_ifft_log.m[index][0] = m_ifft_log.m[index][0];
                m_fold_ifft_log.m[index][1] = 0;
                m_fold_ifft_log.m[index][0] += m_ifft_log.m[m_fft_size - index][0];
            }
            if (index == m_fft_size/2) {
                m_fold_ifft_log.m[index][0] = m_ifft_log.m[index][0];
            }
            if (index > m_fft_size/2) {
                m_fold_ifft_log.m[index][0] = 0;
                m_fold_ifft_log.m[index][1] = 0;
            }
        }

        DBG_COMPLEX_VECTOR("fold:", m_fold_ifft_log.m, m_fft_size)

        // 3 -> 5
        m_dft1.fft (m_fold_ifft_log, m_fft_fold_ifft_log);

        for (size_t index = 0; index < m_fft_size; ++index) {
            std::complex<float> c(m_fft_fold_ifft_log.m[index][0], m_fft_fold_ifft_log.m[index][1]);
            std::complex<float> c2 = std::exp(c);
            m_exp_fft_fold_ifft_log.m[index][0] = std::real(c2);
            m_exp_fft_fold_ifft_log.m[index][1] = std::imag(c2);
        }

        DBG_COMPLEX_VECTOR("exp:", m_exp_fft_fold_ifft_log.m, m_fft_size)

        // 3 -> 4
        m_dft1.ifft (m_exp_fft_fold_ifft_log, m_ifft_exp_fft_fold_ifft_log);

        for (size_t index = 0; index < m_fft_size; ++index)
        {
            m_ifft_exp_fft_fold_ifft_log.m[index][0] /= m_fft_size;
        }

        for (size_t index = 0; index < m_fft_size; ++index)
        {
            m_minimum_phase_response[index] = m_ifft_exp_fft_fold_ifft_log.m[index][0];
        }

        DBG_REAL_VECTOR("linear phase response:", &m_linear_phase_response[0], m_fft_size)
        
        DBG_REAL_VECTOR("minimum phase response:", &m_minimum_phase_response[0], m_fft_size)
    }

protected:
    void add_frames_to_buffer (size_t buffer_index, const float *buffer, size_t number_of_samples)
    {
        dft_buffer &buffer1 = (0 == buffer_index) ? m_buffer11 : m_buffer21;
        dft_buffer &buffer2 = (0 == buffer_index) ? m_buffer12 : m_buffer22;

        dft_buffer &spectrum = (0 == buffer_index) ? m_spectrum1 : m_spectrum2;

        size_t &number_of_ffts = (0 == buffer_index) ? m_number_of_ffts1 : m_number_of_ffts2;

        int &buffer_head1 = 0 == (buffer_index) ? m_buffer_head11 : m_buffer_head21;
        int &buffer_head2 = 0 == (buffer_index) ? m_buffer_head12 : m_buffer_head22;

        for (size_t sample_index = 0; sample_index < number_of_samples; ++sample_index)
        {
            buffer1.m[buffer_head1][0] = buffer[sample_index] * m_window[buffer_head1];
            buffer2.m[buffer_head2][0] = buffer[sample_index] * m_window[buffer_head2];

            ++buffer_head1;
            ++buffer_head2;

            if (buffer_head1 >= (int)m_fft_size) {
                    // DBG("execute " << buffer_index << "-1\n")
                    m_dft1.fft (buffer1, spectrum);
                    for (size_t index = 0; index < m_fft_size; ++index) {
                            spectrum.m[index][0] += sqrt(pow(spectrum.m[index][0], 2) + pow(spectrum.m[index][1], 2));
                            spectrum.m[index][1] = 0;
                    }
                    ++number_of_ffts;
                    buffer_head1 = 0;
            }

            if (buffer_head2 >= (int)m_fft_size) {
                    // DBG("execute " << buffer_index << "-2\n")
                    m_dft1.fft (buffer2, spectrum);
                    for (size_t index = 0; index < m_fft_size; ++index) {
                            spectrum.m[index][0] += sqrt(pow(spectrum.m[index][0], 2) + pow(spectrum.m[index][1], 2));
                            spectrum.m[index][1] = 0;
                    }
                    ++number_of_ffts;
                    buffer_head2 = 0;
            }
        }
    }

    void reset_buffer (size_t buffer_index)
    {
        dft_buffer &buffer1 = (0 == buffer_index) ? m_buffer11 : m_buffer21;
        dft_buffer &buffer2 = (0 == buffer_index) ? m_buffer12 : m_buffer22;

        dft_buffer &spectrum = (0 == buffer_index) ? m_spectrum1 : m_spectrum2;

        size_t &number_of_ffts = (0 == buffer_index) ? m_number_of_ffts1 : m_number_of_ffts2;

        int &buffer_head1 = 0 == (buffer_index) ? m_buffer_head11 : m_buffer_head21;
        int &buffer_head2 = 0 == (buffer_index) ? m_buffer_head12 : m_buffer_head22;

        memset(buffer1.m, 0, sizeof(float) * m_fft_size * 2);
        memset(buffer2.m, 0, sizeof(float) * m_fft_size * 2);
        memset(spectrum.m, 0, sizeof(float) * m_fft_size * 2);

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
