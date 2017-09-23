#include "ba4.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    std::vector< std::string > reads;
    if( argc > 1 )
        reads = rosalind::io::getFileLines( argv[1] );
    else reads = rosalind::io::readInputStream();

    std::cout << rosalind::io::join(  rosalind::ba4::substringEncodingAA(
                                          reads.front().cbegin(),
                                          reads.front().cend() ,
                                          reads.at( 1 )) ,
                                      "\n" );
    return 0;
}
