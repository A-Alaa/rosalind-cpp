#include "ba1.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    auto minmax = rosalind::ba1::skewDiagram( input[0] );

    for ( auto minIdx : minmax.first.first )
        std::cout << minIdx << " ";
    return 0;
}
