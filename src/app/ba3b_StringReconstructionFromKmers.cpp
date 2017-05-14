#include "ba3.hpp"

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    std::cout << rosalind::ba3::reconstructStringFromKmers( input.begin() ,
                                                            input.end())
              << std::endl;
    return 0;
}
