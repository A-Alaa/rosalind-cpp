#include "ba4.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    std::string peptide;
    if( argc > 1 )
        peptide = rosalind::io::getFileLines( argv[1] ).front();
    else peptide = rosalind::io::readInputStream().front();

    std::cout << rosalind::io::join(
                     rosalind::io::asStringsVector(
                         rosalind::ba4::theoriticalLinearSpectrum(
                             peptide.cbegin() , peptide.cend())), " " );
    return 0;
}
