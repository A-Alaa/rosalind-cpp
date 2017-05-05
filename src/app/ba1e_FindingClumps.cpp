#include "Common.hh"
#include <iostream>

int main(int argc, char *argv[])
{
    if( argc == 4 )
    {
        auto input = rosalind::io::readInputStream();
        auto clumps = rosalind::basic::findClumps( input[0] ,
                atoi( argv[1]) , atoi( argv[2]) , atoi( argv[3]));

        for ( auto c : clumps )
            std::cout << c << " ";
        return 0;
    }
    else exit( EXIT_FAILURE );
}
