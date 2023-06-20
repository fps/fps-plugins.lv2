#include "stereo_decorrelation.h"

#include <iostream>

int main (int argc, char *argv[])
{
    if (argc < 5)
    {
        std::cout << "usage: test_white_noise_burst length seed decay whiten\n";
        return EXIT_FAILURE;
    }

    size_t length = atoi (argv[1]);
    int seed = atoi (argv[2]);
    float decay = atof (argv[3]);
    bool whiten = atoi (argv[4]) != 0;

    std::vector<float> data (length, 0);
    create_exponential_white_noise_burst (seed, decay, whiten, data);

    for (size_t index = 0; index < data.size (); ++index)
    {
        std::cout << index << " " << data[index] << "\n";
    }
}
