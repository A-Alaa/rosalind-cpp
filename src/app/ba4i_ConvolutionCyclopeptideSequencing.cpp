#include "ba4.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    std::vector< std::string > input;
    if( argc > 1 )
        input = rosalind::io::getFileLines( argv[1] );
    else input = rosalind::io::readInputStream();
    int m = std::stoi( input.front());
    int n = std::stoi( input.at( 1 ));
    std::vector< int > spectrum;
    for( const std::string &strMass : rosalind::io::split( input.at( 2 ) , " "))
        spectrum.push_back( std::stoi( strMass ));

    std::vector< std::pair< std::vector< uint8_t > , int >>
            sequence =
            rosalind::ba4::convoultionCyclopeptideSequenceing( m , n , spectrum.cbegin() ,
                                                               spectrum.cend());

    std::cout << rosalind::io::join(
                     rosalind::io::asStringsVector( sequence.front().first ) , "-");
    return 0;
}
