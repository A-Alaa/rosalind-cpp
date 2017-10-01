#include "ba4.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    std::vector< std::string > input;
    if( argc > 1 )
        input = rosalind::io::getFileLines( argv[1] );
    else input = rosalind::io::readInputStream();
    std::vector< int > spectrum;
    for( const std::string &strMass : rosalind::io::split( input.at( 1 ) , " "))
        spectrum.push_back( std::stoi( strMass ));
    std::cout << rosalind::ba4::scoreCyclicPeptideSpectrum( input.front() ,
                                                            spectrum.cbegin() ,
                                                            spectrum.cend() )
              << std::endl;
    return 0;
}
