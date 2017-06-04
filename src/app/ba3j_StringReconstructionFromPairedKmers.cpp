#include "ba3.hpp"

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    unsigned int d = atoi( rosalind::io::split( input.front() , " ")[1].c_str() );
    std::cout << rosalind::ba3::reconstructStringFromSparsePairedKmers(
                     input.cbegin() + 1 , input.cend(), d );
    return 0;
}
