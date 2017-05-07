#include "Common.hh"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    std::cout << rosalind::basic::hammingDistance( input[0] , input[1] );
    return 0;
}
