#include "ba3.hpp"

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    for( const auto &overlap : rosalind::ba3::constructDeBruijnGraph( input.cbegin() , input.cend()))
        std::cout << overlap.first << " -> " << rosalind::io::join( overlap.second , "," )
                  << std::endl;
    return 0;
}
