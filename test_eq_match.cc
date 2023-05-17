#define EQ_MATCH_FLOAT float

#include <sndfile.h>
#include <stdexcept>
#include <iostream>

#include "eq_match.h"

#define FFT_SIZE 1024

int main () 
{
  {
    SF_INFO input_info;
    input_info.format = 0;

    SF_INFO output_info;
    output_info.format = 0;

    SNDFILE *input_file = sf_open ("input.wav", SFM_READ, &input_info);
    SNDFILE *output_file = sf_open ("output.wav", SFM_READ, &output_info);

    if (0 == input_file || 0 == output_file) throw std::runtime_error ("Failed to open one or the other file");

    sf_count_t input_count = sf_seek (input_file, 0, SEEK_END);
    sf_count_t output_count = sf_seek (output_file, 0, SEEK_END);

    std::cout << input_count << " " << output_count << "\n";

    std::vector<float> input (input_count);
    std::vector<float> output (output_count);

    sf_seek (input_file, 0, SEEK_SET);
    sf_seek (output_file, 0, SEEK_SET);

    sf_count_t input_read = sf_read_float (input_file, &input[0], input_count);
    sf_count_t output_read = sf_read_float (output_file, &output[0], output_count);

    std::cout << input_read << " " << output_read << "\n";


    eq_match match (FFT_SIZE, input_info.samplerate);
    match.reset_buffer1 ();
    match.reset_buffer2 ();

    for (size_t index = 0; index < (input_count - FFT_SIZE); index += (FFT_SIZE/2))
    {
      match.add_frames_to_buffer1 (&input[index], FFT_SIZE);
      match.add_frames_to_buffer2 (&output[index], FFT_SIZE);
    }

    match.calculate_response ();

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
  }
  // fftwf_cleanup ();
}

