#include "Common.hh"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    for( auto frequency : rosalind::basic::stringFrequencyArray( input ))
        std::cout << frequency << " ";
    return 0;
}
