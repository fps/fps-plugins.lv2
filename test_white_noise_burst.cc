#include "stereo_decorrelation.h"

#include <iostream>

int main (int argc, char *argv[])
{
    std::vector<float> data (1000, 0);
    create_exponential_white_noise_burst (0, 100, false, data);

    std::cout << "non-whitened: ";
    for (size_t index = 0; index < data.size (); ++index)
    {
        std::cout << data[index] << " ";
    }
    std::cout << "\n";

    create_exponential_white_noise_burst (0, 100, true, data);

    std::cout << "whitened: ";
    for (size_t index = 0; index < data.size (); ++index)
    {
        std::cout << data[index] << " ";
    }
    std::cout << "\n";

}
