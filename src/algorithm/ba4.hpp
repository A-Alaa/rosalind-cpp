#ifndef BA4_HPP
#define BA4_HPP
#include "ba3.hpp"

namespace rosalind
{
namespace ba4
{
constexpr std::array< const char , 4 > gacu = { 'G' , 'A' , 'C' , 'U' };
constexpr std::array< const char , 4 > gact = { 'G' , 'A' , 'C' , 'T' };

std::array< int , byteCapacity > codeGACUT([]{
    std::array< int , byteCapacity > codes;
    std::generate_n( codes.begin() , byteCapacity , [](){ return -1;});
    codes['G'] = 0;
    codes['A'] = 1;
    codes['C'] = 2;
    codes['U'] = 3;
    codes['T'] = 3;
    return codes;
}());

static char codonToAcid( std::size_t index )
{
    assert( index < 64 );
    static const std::array< const char, 65 > codonToAcid =
    { "GGGG"
      "EEDD"
      "AAAA"
      "VVVV"
      "RRSS"
      "KKNN"
      "TTTT"
      "MIII"
      "RRRR"
      "QQHH"
      "PPPP"
      "LLLL"
      "W*CC"
      "**YY"
      "SSSS"
      "LLFF" };
    return codonToAcid.at( index );
}

static char codonToAcid( char b1, char b2, char b3 )
{
    return codonToAcid( codeGACUT[ b1 ] * powi< 4 >( 2 ) +
                        codeGACUT[ b2 ] * powi< 4 >( 1 ) +
                        codeGACUT[ b3 ] * powi< 4 >( 0 ));
}


template< typename CharIt >
std::string cdogma( const CharIt first , const CharIt last )
{
    const int L = std::distance( first , last );
    assert(  L > 0 && L % 3 == 0 );
    std::string aa;
    aa.reserve( L / 3 );
    for( auto it = first ; it != last ; it+=3 )
        aa.push_back( codonToAcid( *it ,
                                   *std::next(it) ,
                                   *std::next(it,2)));
    return aa;
}

std::string cdogma( const std::string &s )
{
    return cdogma( s.cbegin() , s.cend());
}

template< typename CharIt >
std::string translateRNA2AA( const CharIt first , const CharIt last )
{
    assert( std::all_of( first , last , []( char c ){
        return std::any_of( ugca.cbegin() , ugca.cend() ,
                            [c]( char b ){ return b == c ;});
    }));
    return cdogma( first , last );
}

template< typename CharIt >
std::vector< std::string >
substringEncodingAA( CharIt genomeFirstIt ,
                     CharIt genomeLastIt ,
                     const std::string &aa )
{
    const size_t k = aa.size() * 3;
    std::vector< std::string > results;
    if( std::distance( genomeFirstIt , genomeLastIt ) >= k )
        for( auto it1 = genomeFirstIt , it2 = std::next( genomeFirstIt , k ) ;;
             ++it1 , ++it2 )
        {
            const std::string kmer{ it1 , it2 }; // TODO: use std++17 std::string_view.
            if( aa == cdogma( kmer ) ||
                    aa == cdogma( ba1::complementSequence( kmer )))
                results.emplace_back( it1 , it2 );
            if( it2 == genomeLastIt )
                break;
        }
    return results;
}

}
}

#endif // BA4_HPP
