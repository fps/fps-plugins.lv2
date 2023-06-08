#include <sndfile.h>
#include <stdexcept>
#include <iostream>
#include <FFTConvolver/FFTConvolver.h>

#include "eq_match.h"

int main (int argc, char *argv[])
{
    if (argc < 7)
    {
        std::cerr << "Usage: test_eq_match FFT_SIZE input_file1 input_file2 linear_phase_response_output_file minimum_phase_response_output_file matched_output_file\n";
        return EXIT_FAILURE;
    }

    const int FFT_SIZE = atoi (argv[1]);

    std::cerr << "Reading input1 files...\n";

    SF_INFO input1_info;
    input1_info.format = 0;

    SF_INFO input2_info;
    input2_info.format = 0;

    SNDFILE *input1_file = sf_open (argv[2], SFM_READ, &input1_info);
    SNDFILE *input2_file = sf_open (argv[3], SFM_READ, &input2_info);

    if (0 == input1_file) throw std::runtime_error ("Failed to open inputfile1");

    if (0 == input2_file) throw std::runtime_error ("Failed to open inputfile2");

    if (input1_info.samplerate != input2_info.samplerate) throw std::runtime_error ("inputfile1 and inputfile2 have different samplerates");

    if (1 != input1_info.channels || 1 != input2_info.channels) throw std::runtime_error ("Channels != 1");

    sf_count_t input1_count = sf_seek (input1_file, 0, SEEK_END);
    sf_count_t input2_count = sf_seek (input2_file, 0, SEEK_END);

    // std::cout << input1_count << " " << input2_count << "\n";

    std::vector<float> input1 (input1_count);
    std::vector<float> input2 (input2_count);

    sf_seek (input1_file, 0, SEEK_SET);
    sf_seek (input2_file, 0, SEEK_SET);

    sf_count_t input1_read = sf_read_float (input1_file, &input1[0], input1_count);
    sf_count_t input2_read = sf_read_float (input2_file, &input2[0], input2_count);

    sf_close (input1_file);
    sf_close (input2_file);

    std::cerr << "Calculating spectra...\n";

    eq_match match (FFT_SIZE);

    match.add_frames_to_buffer1 (&input1[0], input1_count);
    match.add_frames_to_buffer2 (&input2[0], input1_count);

    std::cerr << "Calculating responses...\n";
    match.calculate_response ();

#ifdef DEBUG    
    std::cout << "linear phase response: ";
    for (size_t index = 0; index < FFT_SIZE; ++index)
    {
        std::cout << match.m_linear_phase_response[index] << " ";
    }
    std::cout << "\n";

    std::cout << "minimum phase response: ";
    for (size_t index = 0; index < FFT_SIZE; ++index)
    {
        std::cout << match.m_minimum_phase_response[index] << " ";
    }
    std::cout << "\n";
#endif
    
    std::cerr << "Processing files...\n";

    fftconvolver::FFTConvolver convolver;
    convolver.init (64, &match.m_minimum_phase_response[0], FFT_SIZE);
    // convolver.init (64, &match.m_linear_phase_response[0], FFT_SIZE);
    std::vector<float> minimum_phase_output (input1_count);

    convolver.process (&input1[0], &minimum_phase_output[0], input1_count);

    std::cerr << "Writing output files...\n";

    SF_INFO linear_phase_response_info;
    SF_INFO minimum_phase_response_info;
    SF_INFO matched_output_info;

    linear_phase_response_info.channels = 1;
    minimum_phase_response_info.channels = 1;
    matched_output_info.channels = 1;

    linear_phase_response_info.samplerate = input1_info.samplerate;
    minimum_phase_response_info.samplerate = input1_info.samplerate;
    matched_output_info.samplerate = input1_info.samplerate;

    linear_phase_response_info.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;
    minimum_phase_response_info.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;
    matched_output_info.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;

    SNDFILE *linear_phase_output_file = sf_open (argv[4], SFM_WRITE, &linear_phase_response_info);

    SNDFILE *minimum_phase_output_file = sf_open (argv[5], SFM_WRITE, &minimum_phase_response_info);

    SNDFILE *matched_output_file = sf_open (argv[6], SFM_WRITE, &matched_output_info);

    if (0 == linear_phase_output_file) throw std::runtime_error ("Failed to open linear_phase_response output file: " + std::string (argv[4]));

    if (0 == minimum_phase_output_file) throw std::runtime_error ("Failed to open minimum_phase_response output file");

    if (0 == matched_output_file) throw std::runtime_error ("Failed to open matched output file");

    sf_write_float (linear_phase_output_file, &match.m_linear_phase_response[0], FFT_SIZE);

    sf_write_float (minimum_phase_output_file, &match.m_minimum_phase_response[0], FFT_SIZE);

    sf_write_float (matched_output_file, &minimum_phase_output[0], input1_count);

    sf_close (linear_phase_output_file);
    sf_close (minimum_phase_output_file);
    sf_close (matched_output_file);

    std::cerr << "Done.\n";
}

