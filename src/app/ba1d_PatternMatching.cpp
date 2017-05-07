#include "Common.hh"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    auto occurances = rosalind::basic::patternMatching( input[1] , input[0] );
    for( auto idx : occurances )
        std::cout << idx << " ";
    return 0;
}
