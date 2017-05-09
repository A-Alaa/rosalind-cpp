#include "ba1.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    auto kmers = rosalind::ba1::frequentWordsWithMismatches( input );
    std::cout << rosalind::io::join( kmers , " ");
    return 0;
}
