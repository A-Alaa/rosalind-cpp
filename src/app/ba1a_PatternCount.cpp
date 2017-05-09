#include "ba1.hpp"

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    std::cout << rosalind::ba1::patternCount( input[ 0 ] , input[ 1 ]);
    return 0;
}
