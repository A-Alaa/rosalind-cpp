#include "Common.hh"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    auto mostFrequentKmers =
            rosalind::basic::frequentWordsBruteForce( input );
    std::cout << rosalind::io::join( mostFrequentKmers , " ");
    return 0;
}
