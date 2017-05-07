#include "Common.hh"
#include <iostream>

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    auto indices = rosalind::basic::approximatePatternMatching(
                input[1] , input[0] ,  atoi( input[2].c_str() ));
    for ( auto idx : indices )
        std::cout << idx << " ";
    return 0;
}
