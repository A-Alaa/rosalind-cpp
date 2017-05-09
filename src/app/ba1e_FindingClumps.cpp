#include "ba1.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    if( argc == 4 )
    {
        auto input = rosalind::io::readInputStream();
        auto clumps = rosalind::ba1::findClumps( input );

        for ( auto c : clumps )
            std::cout << c << " ";
        return 0;
    }
    else exit( EXIT_FAILURE );
}
