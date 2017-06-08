#ifndef COMMON_HH
#define COMMON_HH

// STL containers
#include <vector>
#include <list>
#include <array>
#include <numeric>
#include <queue>
#include <string>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <iterator>

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
#include <random>

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

template< size_t Base  >
auto powi( uint16_t exponent )
{
    if( exponent == 0 ) return size_t{1} ;
    return  Base * powi< Base >( exponent - 1 );
}


/**
 * @brief random_element
 * credits: http://stackoverflow.com/a/6943003
 * @param begin
 * @param end
 * @return
 */
template <typename I>
I randomElement( I begin, I end )
{
    using UIntType = unsigned long;
    const UIntType n = std::distance(begin, end);
    const UIntType divisor = (UIntType(RAND_MAX) + 1) / n;

    unsigned long k;
    do { k = std::rand() / divisor; } while (k >= n);

    std::advance(begin, k);
    return begin;
}

template <typename I>
auto randomElementSampler( I begin, I end )
{
    static std::mt19937 rng( std::random_device{}());
    const auto n = std::distance( begin , end );
    return [n,begin]() -> I
    {
        std::uniform_int_distribution< decltype(n) > dist{ 0 , n - 1 };
        return std::next( begin , dist( rng ));
    };
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

auto split( const std::string &s , std::string delim  )
{
    std::vector< std::string > tokens;
    size_t last = 0; size_t next = 0;
    while ((next = s.find(delim, last )) != std::string::npos)
    {
        tokens.push_back( s.substr(last , next - last));
        last = next + 1;
    }
    last += delim.length() - 1;
    if( last < s.length() )
        tokens.push_back( s.substr( last , std::string::npos ));
    return tokens;
}

template< typename SeqIt >
std::string join( SeqIt first , SeqIt last , const std::string &sep )
{
    auto binaryJoinString = [sep]( const std::string &a , const std::string &b ) -> std::string
    {
        return  a + ((a.length() > 0) ? sep : "") +  b ;
    };
    return std::accumulate( first , last ,
                            std::string() , binaryJoinString  );
}

template< typename Container = std::vector< std::string >>
std::string join( const Container &container ,
                  const std::string &sep )
{
    return join( container.cbegin() , container.cend() , sep );
}



auto
findWithMismatches( const std::string &str ,
                    const std::string &substr,
                    unsigned int d )
{
    assert( substr.size() <= str.size() && substr.size() > d );
    auto cStr = str.data();
    auto cSubstr = substr.data();
    auto k = substr.size();

    auto occurance = [k,d,&cStr,&cSubstr]( std::string::size_type index )
    {
        decltype(d) mismatches = 0;
        decltype(index) i = 0;
        while( mismatches <= d && i < k )
            mismatches += cStr[ index + i ] != cSubstr[ i++ ];
        return mismatches <= d;
    };

    for( std::string::size_type i = 0, until = str.size() - substr.size() + 1;
         i <  until ; i++ )
        if( occurance( i ))
            return i;

    return std::string::npos;
}
}
}
#endif // COMMON_HH
