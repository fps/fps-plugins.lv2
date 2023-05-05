#ifndef FPS_PLUGINS_EQ_MATCH_HH
#define FPS_PLUGINS_EQ_MATCH_HH

/**
 * This whole class is not realtime safe!
 *
 * How to:
 * (1) Construct an object of type eq_match
 * (2) Call add_frames_to_buffer1 consecutively to estimate spectrum1
 * (3) Call add_frames_to_buffer2 consecutively to estimate spectrum1
 * (4) Call calculate_response ()
 * (5) Get the calculated IRs from m_response (linear phase) or
 *     m_minimum_phase_response (which is minimum phase)
 *
 * You need to provide a macro definition for the floating
 * point type used. e.g.
 *
 *   #define FPS_PLUGINS_EQ_MATCH_FLOATING_POINT_TYPE float
 *
 * or
 *
 *   #define FPS_PLUGINS_EQ_MATCH_FLOATING_POINT_TYPE double
 */

#if FPS_PLUGINS_EQ_MATCH_FLOATING_POINT_TYPE == float
  #define FPS_FFTW_COMPLEX fftwf_complex
  #define FPS_FFTW_PLAN fftwf_plan
  #define FPS_FFTW_PLAN_DFT fftwf_plan_dft_1d
  #define FPS_FFTW_ALLOC_COMPLEX fftwf_alloc_complex
  #define FPS_FFTW_FREE fftwf_free
  #define FPS_FFTW_DESTROY_PLAN fftwf_destroy_plan
  #define FPS_FFTW_EXECUTE fftwf_execute
#elif FPS_PLUGINS_EQ_MATCH_FLOATING_POINT_TYPE == double
  #define FPS_FFTW_COMPLEX fftw_complex
  #define FPS_FFTW_PLAN fftw_plan
  #define FPS_FFTW_PLAN_DFT fftw_plan_dft_1d
  #define FPS_FFTW_ALLOC_COMPLEX fftw_alloc_complex
  #define FPS_FFTW_FREE fftw_free
  #define FPS_FFTW_DESTROY_PLAN fftw_destroy_plan
  #define FPS_FFTW_EXECUTE FPS_FFTW_EXECUTE
#else
  #error "Unexpected value of FPS_PLUGINS_EQ_MATCH_FLOATING_POINT_TYPE"
#endif

#define FPS_FLOAT FPS_PLUGINS_EQ_MATCH_FLOATING_POINT_TYPE

struct eq_match
{
  size_t m_fft_size;
  FPS_FLOAT m_sample_rate;

  std::vector<FPS_FLOAT> m_window;

  // Buffers to accumulate samples
  FPS_FFTW_COMPLEX *m_buffer11;
  int m_buffer_head11;
  FPS_FFTW_COMPLEX *m_buffer12;
  int m_buffer_head12;

  FPS_FFTW_COMPLEX *m_buffer21;
  int m_buffer_head21;
  FPS_FFTW_COMPLEX *m_buffer22;
  int m_buffer_head22;

  FPS_FFTW_COMPLEX *m_spectrum1;
  FPS_FFTW_COMPLEX *m_spectrum2;

  FPS_FFTW_PLAN m_fft_plan11;
  FPS_FFTW_PLAN m_fft_plan12;
  FPS_FFTW_PLAN m_fft_plan21;
  FPS_FFTW_PLAN m_fft_plan22;

  FPS_FFTW_PLAN m_fft_plan;
  FPS_FFTW_PLAN m_ifft_plan;

  FPS_FFTW_COMPLEX *m_fft_buffer1;
  size_t m_number_of_ffts1;

  FPS_FFTW_COMPLEX *m_fft_buffer2;
  size_t m_number_of_ffts2;

  FPS_FFTW_COMPLEX *m_fft_buffer3;
  FPS_FFTW_COMPLEX *m_fft_buffer4;
  FPS_FFTW_COMPLEX *m_fft_buffer5;

  FPS_FFTW_COMPLEX *m_response;
  FPS_FFTW_COMPLEX *m_minimum_phase_response;

  eq_match (size_t fft_size, FPS_FLOAT sample_rate) :
    m_fft_size (fft_size),
    m_sample_rate (sample_rate),
    m_window(fft_size, 0.f)
  {
    m_buffer11 = (FPS_FFTW_COMPLEX*)FPS_FFTW_ALLOC_COMPLEX (fft_size);
    m_buffer12 = (FPS_FFTW_COMPLEX*)FPS_FFTW_ALLOC_COMPLEX (fft_size);

    m_buffer21 = (FPS_FFTW_COMPLEX*)FPS_FFTW_ALLOC_COMPLEX (fft_size);
    m_buffer22 = (FPS_FFTW_COMPLEX*)FPS_FFTW_ALLOC_COMPLEX (fft_size);

    m_buffer_head11 = 0;
    m_buffer_head12 = fft_size/2;
    m_buffer_head21 = 0;
    m_buffer_head22 = fft_size/2;

    m_spectrum1 = (FPS_FFTW_COMPLEX*)FPS_FFTW_ALLOC_COMPLEX (fft_size);
    m_spectrum2 = (FPS_FFTW_COMPLEX*)FPS_FFTW_ALLOC_COMPLEX (fft_size);

    m_fft_buffer1 = (FPS_FFTW_COMPLEX*)FPS_FFTW_ALLOC_COMPLEX (fft_size);
    m_number_of_ffts1 = 0;

    m_fft_buffer2 = (FPS_FFTW_COMPLEX*)FPS_FFTW_ALLOC_COMPLEX (fft_size);
    m_number_of_ffts2 = 0;

    m_fft_buffer3 = (FPS_FFTW_COMPLEX*)FPS_FFTW_ALLOC_COMPLEX (fft_size);
    m_fft_buffer4 = (FPS_FFTW_COMPLEX*)FPS_FFTW_ALLOC_COMPLEX (fft_size);
    m_fft_buffer5 = (FPS_FFTW_COMPLEX*)FPS_FFTW_ALLOC_COMPLEX (fft_size);

    // hann window generation
    for (size_t index = 0; index < fft_size; ++index) {
      m_window[index] = pow(sin(M_PI * index / fft_size), 2);
    }

    m_fft_plan11 = FPS_FFTW_PLAN_DFT(fft_size, &m_buffer11[0], m_fft_buffer1, FFTW_FORWARD, FFTW_ESTIMATE);
    m_fft_plan12 = FPS_FFTW_PLAN_DFT(fft_size, &m_buffer12[0], m_fft_buffer1, FFTW_FORWARD, FFTW_ESTIMATE);

    m_fft_plan21 = FPS_FFTW_PLAN_DFT(fft_size, &m_buffer21[0], m_fft_buffer2, FFTW_FORWARD, FFTW_ESTIMATE);
    m_fft_plan22 = FPS_FFTW_PLAN_DFT(fft_size, &m_buffer22[0], m_fft_buffer2, FFTW_FORWARD, FFTW_ESTIMATE);

    m_fft_plan = FPS_FFTW_PLAN_DFT(fft_size, m_fft_buffer3, m_fft_buffer5, FFTW_FORWARD, FFTW_ESTIMATE);

    m_ifft_plan = FPS_FFTW_PLAN_DFT(fft_size, m_fft_buffer3, m_fft_buffer4, FFTW_BACKWARD, FFTW_ESTIMATE);

    m_response = (FPS_FFTW_COMPLEX*)FPS_FFTW_ALLOC_COMPLEX (fft_size);

    m_minimum_phase_response = (FPS_FFTW_COMPLEX*)FPS_FFTW_ALLOC_COMPLEX (fft_size);
}

  ~eq_match ()
  {
    FPS_FFTW_DESTROY_PLAN (m_fft_plan11);
    FPS_FFTW_DESTROY_PLAN (m_fft_plan12);
    FPS_FFTW_DESTROY_PLAN (m_fft_plan21);
    FPS_FFTW_DESTROY_PLAN (m_fft_plan22);

    FPS_FFTW_DESTROY_PLAN (m_fft_plan);
    FPS_FFTW_DESTROY_PLAN (m_ifft_plan);

    FPS_FFTW_FREE (m_spectrum1);
    FPS_FFTW_FREE (m_spectrum2);

    FPS_FFTW_FREE (m_buffer11);
    FPS_FFTW_FREE (m_buffer12);
    FPS_FFTW_FREE (m_buffer22);
    FPS_FFTW_FREE (m_buffer21);

    FPS_FFTW_FREE (m_fft_buffer1);
    FPS_FFTW_FREE (m_fft_buffer2);
    FPS_FFTW_FREE (m_fft_buffer3);
    FPS_FFTW_FREE (m_fft_buffer4);
    FPS_FFTW_FREE (m_fft_buffer5);

    FPS_FFTW_FREE (m_response);
    FPS_FFTW_FREE (m_minimum_phase_response);
  }


  void add_frames_to_buffer1 (const FPS_FLOAT *buffer, size_t number_of_samples)
  {
    add_frames_to_buffer (0, buffer, number_of_samples);
  }

  void add_frames_to_buffer2 (const FPS_FLOAT *buffer, size_t number_of_samples)
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

  void calculate_response ()
  {
    std::cout << "spectrum1: ";
    for (size_t index = 0; index < m_fft_size; ++index) std::cout << m_spectrum1[index][0] << " ";
    std::cout << "\nspectrum2: ";
    for (size_t index = 0; index < m_fft_size; ++index) std::cout << m_spectrum2[index][0] << " ";
    std::cout << "\n\n";

    for (size_t index = 0; index < m_fft_size; ++index) {
      m_fft_buffer3[index][0] = ((m_spectrum2[index][0] / m_number_of_ffts2) / (m_spectrum1[index][0] / m_number_of_ffts1));
      m_fft_buffer3[index][1] = 0;
    }

    std::cout << "ratio: ";
    for (size_t index = 0; index < m_fft_size; ++index) std::cout << m_fft_buffer3[index][0] << " ";
    std::cout << "\n\n";

    // exp(fft(fold(ifft(log(s)))))
    for (size_t index = 0; index < m_fft_size; ++index) {
      m_fft_buffer3[index][0] = log(m_fft_buffer3[index][0]);
    }

    std::cout << "log: ";
    for (size_t index = 0; index < m_fft_size; ++index) std::cout << m_fft_buffer3[index][0] << " ";
    std::cout << "\n\n";

    std::cout << m_fft_buffer3[0][0] << "\n";

    // 3 -> 4
    FPS_FFTW_EXECUTE(m_ifft_plan);

    for (size_t index = 0; index < m_fft_size; ++index)
    {
      m_fft_buffer4[index][0] /= m_fft_size;
      m_response[index][0] = m_fft_buffer4[index][0];
    }

    std::cout << "ifft: ";
    for (size_t index = 0; index < m_fft_size; ++index) std::cout << m_fft_buffer4[index][0] << " ";
    std::cout << "\n\n";

    std::cout << m_fft_buffer4[0][0] << " " << m_fft_buffer4[0][1] << "\n";

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

    std::cout << "fold: ";
    for (size_t index = 0; index < m_fft_size; ++index) std::cout << m_fft_buffer3[index][0] << " ";
    std::cout << "\n\n";

    std::cout << m_fft_buffer3[0][0] << "\n";

    // 3 -> 5
    FPS_FFTW_EXECUTE (m_fft_plan);

    std::cout << m_fft_buffer5[0][0] << " " << m_fft_buffer5[0][1] << "\n";

    for (size_t index = 0; index < m_fft_size; ++index) {
      std::complex<double> c(m_fft_buffer5[index][0], m_fft_buffer5[index][1]);
      std::complex<double> c2 = std::exp(c);
      m_fft_buffer3[index][0] = std::real(c2);
      m_fft_buffer3[index][1] = std::imag(c2);
    }

    std::cout << "exp: ";
    for (size_t index = 0; index < m_fft_size; ++index) std::cout << m_fft_buffer3[index][0] << " ";
    std::cout << "\n\n";
    for (size_t index = 0; index < m_fft_size; ++index) std::cout << m_fft_buffer3[index][1] << " ";
    std::cout << "\n\n";

    std::cout << m_fft_buffer3[0][0] << "\n";

    // 3 -> 4
    FPS_FFTW_EXECUTE (m_ifft_plan);
    for (size_t index = 0; index < m_fft_size; ++index) m_fft_buffer4[index][0] /= m_fft_size;

    // ifftshift
    // for (size_t index = 0; index < m_fft_size; ++index) response[index][0] = fft_buffer4[(index + m_fft_size/2) % m_fft_size][0];
    for (size_t index = 0; index < m_fft_size; ++index) m_minimum_phase_response[index][0] = m_fft_buffer4[index][0];

    std::cout << "response: ";
    for (size_t index = 0; index < m_fft_size; ++index) std::cout << m_minimum_phase_response[index][0] << " ";
    std::cout << "\n";
    for (size_t index = 0; index < m_fft_size; ++index) std::cout << m_minimum_phase_response[index][1] << " ";
    std::cout << "\n";
    // for (size_t index = 0; index < m_fft_size/2; ++index) response[index + m_fft_size/2][0] = 0;
  }

protected:
  void add_frames_to_buffer (size_t buffer_index, const FPS_FLOAT *buffer, size_t number_of_samples)
  {
    FPS_FFTW_COMPLEX *buffer1 = (0 == buffer_index) ? m_buffer11 : m_buffer21;
    FPS_FFTW_COMPLEX *buffer2 = (0 == buffer_index) ? m_buffer12 : m_buffer22;

    FPS_FFTW_PLAN fft_plan1 = (0 == buffer_index) ? m_fft_plan11 : m_fft_plan21;
    FPS_FFTW_PLAN fft_plan2 = (0 == buffer_index) ? m_fft_plan12 : m_fft_plan22;

    FPS_FFTW_COMPLEX *fft_buffer = (0 == buffer_index) ? m_fft_buffer1 : m_fft_buffer2;

    FPS_FFTW_COMPLEX *spectrum = (0 == buffer_index) ? m_spectrum1 : m_spectrum2;

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
          // std::cout << "execute " << buffer_index << "-1\n";
          FPS_FFTW_EXECUTE (fft_plan1);
          for (size_t index = 0; index < m_fft_size; ++index) {
              spectrum[index][0] += sqrt(pow(fft_buffer[index][0], 2) + pow(fft_buffer[index][1], 2));
              spectrum[index][1] = 0;
          }
          ++number_of_ffts;
          buffer_head1 = 0;
      }

      if (buffer_head2 >= (int)m_fft_size) {
          // std::cout << "execute " << buffer_index << "-2\n";
          FPS_FFTW_EXECUTE (fft_plan2);
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
    FPS_FFTW_COMPLEX *buffer1 = (0 == buffer_index) ? m_buffer11 : m_buffer21;
    FPS_FFTW_COMPLEX *buffer2 = (0 == buffer_index) ? m_buffer12 : m_buffer22;

    FPS_FFTW_COMPLEX *spectrum = (0 == buffer_index) ? m_spectrum1 : m_spectrum2;

    size_t &number_of_ffts = (0 == buffer_index) ? m_number_of_ffts1 : m_number_of_ffts2;

    int &buffer_head1 = 0 == (buffer_index) ? m_buffer_head11 : m_buffer_head21;
    int &buffer_head2 = 0 == (buffer_index) ? m_buffer_head12 : m_buffer_head22;

    memset(buffer1, 0, sizeof(FPS_FLOAT) * m_fft_size * 2);
    memset(buffer2, 0, sizeof(FPS_FLOAT) * m_fft_size * 2);
    memset(spectrum, 0, sizeof(FPS_FLOAT) * m_fft_size * 2);
    buffer_head1 = 0;
    buffer_head2 = m_fft_size/2;
    number_of_ffts = 0;
  }


};

#undef FPS_FFTW_COMPLEX
#undef FPS_FFTW_PLAN
#undef FPS_FFTW_PLAN_DFT
#undef FPS_FFTW_ALLOC_COMPLEX
#undef FPS_FFTW_FREE
#undef FPS_FFTW_DESTROY_PLAN
#undef FPS_FFTW_EXECUTE


#endif
