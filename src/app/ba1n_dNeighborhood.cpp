#include "ba1.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    for( const auto &mutant : rosalind::ba1::dNeighborhood( input ))
        std::cout << mutant << std::endl;
    return 0;
}
