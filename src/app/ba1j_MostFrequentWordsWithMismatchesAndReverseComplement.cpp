#include "Common.hh"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    auto kmers = rosalind::basic::frequentWordsWithMismatchesAndReverseComplement( input );
    std::cout << rosalind::io::join( kmers , " ");
    return 0;
}
