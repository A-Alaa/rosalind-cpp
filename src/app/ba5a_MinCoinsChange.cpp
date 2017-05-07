#include "Common.hh"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    std::cout << rosalind::basic::minimumCoinsChange( input );
    return 0;
}
