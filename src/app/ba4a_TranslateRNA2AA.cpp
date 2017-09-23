#include "ba4.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    std::vector< std::string > reads;
    if( argc > 1 )
        reads = rosalind::io::getFileLines( argv[1] );
    else reads = rosalind::io::readInputStream();

    std::cout << rosalind::ba4::translateRNA2AA( reads.front().cbegin(),
                                                     reads.front().cend());
    return 0;
}
