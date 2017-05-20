#include "ba3.hpp"

int main(int argc, char *argv[])
{
    auto input = rosalind::io::readInputStream();
    std::vector< std::string > out;
    auto _ = rosalind::ba3::findEulerianCycle( input.begin() , input.end());
    std::transform( std::begin( _ ) , std::end( _ ) ,
                    std::inserter( out , std::end( out )) ,
                    []( const rosalind::ba3::Graph< int >::Vertex &v )
    {
        return std::to_string( v.data());
    });

    std::cout << rosalind::io::join(  out , "->")
              << std::endl;
    return 0;
}
