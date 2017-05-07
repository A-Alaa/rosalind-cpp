#include "Common.hh"

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    std::cout << rosalind::basic::patternCount( input[ 0 ] , input[ 1 ]);
    return 0;
}
