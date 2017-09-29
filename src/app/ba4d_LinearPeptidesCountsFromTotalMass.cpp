#include "ba4.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    std::string mass;
    if( argc > 1 )
        mass = rosalind::io::getFileLines( argv[1] ).front();
    else mass = rosalind::io::readInputStream().front();

    std::cout << rosalind::ba4::spectrumsCountFromTotalMass(
                     std::stoi( mass ));
    return 0;
}
