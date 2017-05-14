#include "ba3.hpp"

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    for( const auto &overlap : rosalind::ba3::overlapStringsGraph( input.cbegin() , input.cend()))
        std::cout << overlap.first << " -> " << overlap.second
                  << std::endl;
    return 0;
}
