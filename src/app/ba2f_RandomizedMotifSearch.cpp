#include "ba2.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    for( const auto &motif : rosalind::ba2::randomizedMotifSearch( input ))
        std::cout << motif << std::endl;
    return 0;
}
