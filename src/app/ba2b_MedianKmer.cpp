#include "ba2.hpp"

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    std::cout << rosalind::ba2::medianKmer( input ) << std::endl;
    return 0;
}
