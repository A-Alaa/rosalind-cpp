#include "ba1.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    auto occurances = rosalind::ba1::patternMatching( input );
    for( auto idx : occurances )
        std::cout << idx << " ";
    return 0;
}
