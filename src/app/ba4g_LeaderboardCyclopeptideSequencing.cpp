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
    auto leadingPeptide =
            rosalind::ba4::leaderboardCyclopeptideSequencing(
                std::stoi( input.front()) ,
                spectrum.cbegin() ,
                spectrum.cend()).front();

    std::cout << rosalind::io::join(
                     rosalind::io::asStringsVector( leadingPeptide.first ) , "-");
    return 0;
}
