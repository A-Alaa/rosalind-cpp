#include "ba1.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    for( auto frequency : rosalind::ba1::stringFrequencyArray( input ))
        std::cout << frequency << " ";
    return 0;
}
