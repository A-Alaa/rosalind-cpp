#include "ba3.hpp"

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    for( const auto &kmer : rosalind::ba1::extractKmers( input[1] , atoi( input[0].c_str())))
        std::cout << kmer << std::endl;
    return 0;
}
