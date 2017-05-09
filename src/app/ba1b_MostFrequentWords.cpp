#include "ba1.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    auto mostFrequentKmers =
            rosalind::ba1::frequentWordsBruteForce( input );
    std::cout << rosalind::io::join( mostFrequentKmers , " ");
    return 0;
}
