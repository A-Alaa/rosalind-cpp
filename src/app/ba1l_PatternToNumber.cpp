#include "ba1.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    std::cout << rosalind::ba1::encode< uint64_t >( input[0] );
    return 0;
}
