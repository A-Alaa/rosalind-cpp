#include "ba3.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    std::vector< std::string > reads;
    if( argc > 1 )
        reads = rosalind::io::getFileLines( argv[1] );
    else reads = rosalind::io::readInputStream();

    int k = std::stoi( rosalind::io::split( reads.front() , " ").front());
    int d = std::stoi( rosalind::io::split( reads.front() , " ").at( 1 ));

    std::cout << rosalind::ba3::constructPairedStringsFromGappedPath(
                     reads.cbegin() + 1 , reads.cend() , k , d );
    return 0;
}
