#ifndef BA4_HPP
#define BA4_HPP
#include "ba3.hpp"

namespace rosalind
{
namespace ba4
{
constexpr std::array< const char , 4 > gacu = { 'G' , 'A' , 'C' , 'U' };
constexpr std::array< const char , 4 > gact = { 'G' , 'A' , 'C' , 'T' };

const std::array< int , byteCapacity > codeGACUT([]{
    std::array< int , byteCapacity > codes;
    std::fill( codes.begin() , codes.end() , -1 );
    codes['G'] = 0;
    codes['A'] = 1;
    codes['C'] = 2;
    codes['U'] = 3;
    codes['T'] = 3;
    return codes;
}());

const std::array< uint8_t , byteCapacity > massTable([]{
    std::array< uint8_t , byteCapacity > masses;
    std::fill( masses.begin() , masses.end() , 0 );
    masses['G']=57; masses['V']=99 ; masses['L']=113; masses['E']=129; masses['R']=156;
    masses['A']=71; masses['T']=101; masses['N']=114; masses['M']=131; masses['Y']=163;
    masses['S']=87; masses['C']=103; masses['D']=115; masses['H']=137; masses['W']=186;
    masses['P']=97; masses['I']=113; masses['Q']=128; masses['F']=147; masses['K']=128;
    return masses;
}());

const std::array< uint8_t , 20 > aaMasses([]{
    std::array< uint8_t , 20 > _aaMasses;
    int i = 0;
    std::for_each( massTable.cbegin() , massTable.cend() ,
                   [&_aaMasses,&i]( uint8_t m ){
        if( m > 0 )
            _aaMasses[ i++ ] = m;
    });
    std::sort( _aaMasses.begin() , _aaMasses.end());
    return _aaMasses;
}());

const std::set< uint8_t > aaMassesUnique{ aaMasses.cbegin() , aaMasses.cend()};

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
        return std::any_of( gacu.cbegin() , gacu.cend() ,
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

template< typename CharIt >
std::vector< int >
theoriticalCycloSpectrum( const CharIt aaFirstIt ,
                          const CharIt aaLastIt )
{
    const size_t n = std::distance( aaFirstIt , aaLastIt );
    assert( n > 0 );
    std::vector< int > spectrum;
    spectrum.reserve( n *  n - n + 1 );
    spectrum.push_back( 0 );
    for( auto fragmentSize = 1 ; fragmentSize < n ; ++fragmentSize )
        for( auto seedIt = aaFirstIt; seedIt != aaLastIt ; ++seedIt )
        {
            auto currentIt = seedIt;
            int mass = 0;
            for( auto i = 0 ; i < fragmentSize ; ++i )
            {
                mass += massTable.at( *currentIt );
                if( ++currentIt == aaLastIt )
                    currentIt = aaFirstIt;
            }
            spectrum.push_back( mass );
        }
    spectrum.push_back( std::accumulate( aaFirstIt , aaLastIt , 0 ,
                                         []( int mass , char aa ){
        return mass + massTable.at( aa );
    }));
    return spectrum;
}

uint64_t spectrumsCountFromTotalMass( int mass )
{
    std::vector< uint64_t > peptidesCount( mass + 1 , 0 );
    for( int i = 0 ; i < mass + 1 ; ++i )
        for( const auto aaMass : aaMassesUnique )
        {
            const int parentMass = i - aaMass;
            if( parentMass < 0 )
                break;
            if( parentMass == 0 )
                ++peptidesCount[ i ];
            else
                peptidesCount[ i ] += peptidesCount.at( parentMass );
        }
    return peptidesCount.at( mass );
}

}
}

#endif // BA4_HPP
