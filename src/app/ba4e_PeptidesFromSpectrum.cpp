#include "ba4.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    std::string masses;
    if( argc > 1 )
        masses = rosalind::io::getFileLines( argv[1] ).front();
    else masses = rosalind::io::readInputStream().front();

    const std::vector< std::vector< uint8_t >> numericPeptides =
            rosalind::ba4::massRepresentedPeptidesFromSpectrum( masses );
    std::vector< std::string > out;
    for( const std::vector< uint8_t > &peptide : numericPeptides )
        out.emplace_back(
                    rosalind::io::join(
                        rosalind::io::asStringsVector( peptide ) , "-"));

    std::cout << rosalind::io::join( out , " ") << std::endl;
    return 0;
}
