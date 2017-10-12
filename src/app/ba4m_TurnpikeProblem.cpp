#include "ba4.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    std::vector< std::string > input;
    if( argc > 1 )
        input = rosalind::io::getFileLines( argv[1] );
    else input = rosalind::io::readInputStream();

    std::vector< int > integers;
    for( const std::string &strInteger : rosalind::io::split( input.front() , " "))
        integers.push_back( std::stoi( strInteger ));

    auto reconstructDeltas =
            rosalind::ba4::turnpikeProblem( integers.cbegin() , integers.cend());


    return 0;
}
