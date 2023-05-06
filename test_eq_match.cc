#define FPS_PLUGINS_EQ_MATCH_FLOATING_POINT_TYPE float
#include "eq_match.h"

int main () 
{
  {
    eq_match match (1024, 48000);
    match.reset_buffer1 ();
    match.reset_buffer2 ();

    std::vector<float> buffer (100);

    for (size_t index = 0; index < 1000; ++index)
    {
      match.add_frames_to_buffer1 (&buffer[0], 100);
      match.add_frames_to_buffer2 (&buffer[0], 100);
    }

    match.calculate_response ();
  }
  fftwf_cleanup ();
}

