#include "ba3.hpp"

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    std::cout << rosalind::ba3::reconstructStringFromSparseKmers(
                     input.begin() + 1 , input.end());
    return 0;
}
