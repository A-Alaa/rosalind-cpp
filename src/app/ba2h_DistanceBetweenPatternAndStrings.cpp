#include "ba2.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    std::cout << rosalind::ba2::minimalHammingDistance( input );
    return 0;
}
