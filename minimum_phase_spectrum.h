#ifndef EQ_MATCH_MINIMUM_PHASE_SPECTRUM__
#define EQ_MATCH_MINIMUM_PHASE_SPECTRUM__

#include <cassert>
#include <complex>

#include "dft.h"

void create_minimum_phase_spectrum (const dft_buffer &spectrum_in, dft_buffer &spectrum_out)
{
    assert (spectrum_in.m_size == spectrum_out.m_size);

    const size_t fft_size = spectrum_in.m_size;

    dft the_dft (fft_size);

    // exp(fft(fold(ifft(log(s)))))

    dft_buffer the_log (fft_size);
    for (size_t index = 0; index < fft_size; ++index) {
        std::complex c (spectrum_in.m[index][0], spectrum_in.m[index][1]);
        std::complex c2 = std::log (c);
        the_log.m[index][0] = std::real(c2);
        the_log.m[index][1] = std::imag(c2);
    }

    DBG_COMPLEX_VECTOR("log:", the_log.m, fft_size)

    // 3 -> 4
    dft_buffer ifft_log (fft_size);
    the_dft.ifft (the_log, ifft_log);

    DBG_COMPLEX_VECTOR("ifft:", ifft_log.m, fft_size)

    dft_buffer folded (fft_size);
    // fold
    for (size_t index = 0; index < fft_size; ++index) {
        if (index == 0) {
            folded.m[index][0] = ifft_log.m[index][0];
            folded.m[index][1] = 0;
        }
        if (index > 0 && index < fft_size / 2) {
            folded.m[index][0] = ifft_log.m[index][0];
            folded.m[index][1] = 0;
            folded.m[index][0] += ifft_log.m[fft_size- index][0];
        }
        if (index == fft_size / 2) {
            folded.m[index][0] = ifft_log.m[index][0];
            folded.m[index][1] = 0;
        }
        if (index > fft_size / 2) {
            folded.m[index][0] = 0;
            folded.m[index][1] = 0;
        }
    }

    DBG_COMPLEX_VECTOR("fold:", folded.m, fft_size)

    // 3 -> 5
    dft_buffer fft_folded (fft_size);
    the_dft.fft (folded, fft_folded);

    for (size_t index = 0; index < fft_size; ++index) {
        std::complex<float> c(fft_folded.m[index][0], fft_folded.m[index][1]);
        std::complex<float> c2 = std::exp(c);
        spectrum_out.m[index][0] = std::real(c2);
        spectrum_out.m[index][1] = std::imag(c2);
    }

    DBG_COMPLEX_VECTOR("exp:", spectrum_out.m, fft_size)
}

void create_minimum_phase_response (const std::vector<float> &in, std::vector<float> out)
{

}

#endif
