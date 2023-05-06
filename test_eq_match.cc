#define FPS_PLUGINS_EQ_MATCH_FLOATING_POINT_TYPE float
#include "eq_match.h"

int main () 
{
  {
    eq_match match (1024, 48000);
  }
  fftwf_cleanup ();
}

