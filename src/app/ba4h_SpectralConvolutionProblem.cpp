#include "ba4.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    std::string input;
    if( argc > 1 )
        input = rosalind::io::getFileLines( argv[1] ).front();
    else input = rosalind::io::readInputStream().front();
    std::vector< int > spectrum;
    for( const std::string &strMass : rosalind::io::split( input , " "))
        spectrum.push_back( std::stoi( strMass ));

    auto spectralConvolution =
            rosalind::ba4::spectralConvoultion( spectrum.begin() , spectrum.end());

    std::cout << rosalind::io::join(
                     rosalind::io::asStringsVector( spectralConvolution ) , " ");
    return 0;
}
