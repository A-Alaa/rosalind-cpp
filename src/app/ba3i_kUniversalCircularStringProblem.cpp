#include "ba3.hpp"

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    unsigned int k = std::stoi( input.front());
    std::cout << rosalind::ba3::kUniversalCircularString( k );
    return 0;
}
