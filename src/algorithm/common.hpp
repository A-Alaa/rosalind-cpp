#ifndef COMMON_HH
#define COMMON_HH

// STL containers
#include <vector>
#include <list>
#include <array>
#include <numeric>
#include <queue>
#include <string>
#include <unordered_set>
#include <unordered_map>

// STL streaming
#include <iostream>
#include <sstream>
#include <fstream>

// STL misc
#include <algorithm>
#include <cmath>
#include <typeinfo>
#include <functional>
#include <cassert>

// Qt
//#include <QString>
//#include <QStringRef>

#ifdef _WIN32
#    ifdef LIBRARY_EXPORTS
#        define LIBRARY_API __declspec(dllexport)
#    else
#        define LIBRARY_API __declspec(dllimport)
#    endif
#elif
#    define LIBRARY_API
#endif

#define EXTERN_DLL_EXPORT extern "C" __declspec(dllexport)

namespace rosalind
{

const std::size_t byteCapacity = std::numeric_limits< char >::max() - std::numeric_limits< char >::min();
const std::array< char , 4 > acgt = { 'A' , 'C' , 'G' , 'T' };
const std::array< char , 4 > tgca = { 'T' , 'G' , 'C' , 'A' };
const std::array< std::string , 4 > bpMutants = { "CGT" , "AGT" , "ACT" , "ACG"};

const std::array< int , byteCapacity > codeACGT([]{
    std::array< int , byteCapacity > codes;
    std::generate_n( codes.begin() , byteCapacity , [](){ return -1;});
    codes['A'] = 0;
    codes['C'] = 1;
    codes['G'] = 2;
    codes['T'] = 3;
    return codes;
}());

using IndexType = std::size_t;
using CodeType = std::size_t;
using CountType = std::size_t;
using RosalindIOType = std::vector< std::string >;

template< typename T >
auto powi( T base , uint16_t exponent )
{
    static_assert( std::numeric_limits< T >::is_integer ,
                   "exponent must be of integer type.");
    if( exponent == 0 ) return T( 1 );
    return T( base ) * powi( base , exponent - 1 );
}

namespace io
{

std::vector< std::string >
readInputStream()
{
    std::string line;
    std::vector< std::string > lines;
    while( std::getline( std::cin , line ))
        lines.push_back( line );

    return lines;
}

auto split( const std::string &s , char delim  )
{
    std::stringstream ss( s );
    std::vector< std::string > tokens;
    std::string token;
    while( std::getline( ss , token , delim ))
        tokens.push_back( token );
    return tokens;
}

template< typename Container = std::vector< std::string >>
std::string join( const Container &container ,
                  const std::string &sep )
{
    auto binaryJoinString = [sep]( const std::string &a , const std::string &b ) -> std::string
    {
        return  a + ((a.length() > 0) ? sep : "") +  b ;
    };
    return std::accumulate( container.begin() , container.end() ,
                            std::string() , binaryJoinString  );
}

}
}
#endif // COMMON_HH
