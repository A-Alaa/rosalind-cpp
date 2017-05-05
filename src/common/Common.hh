#ifndef COMMON_HH
#define COMMON_HH

#include <vector>
#include <list>
#include <map>
#include <array>
#include <numeric>
#include <queue>
#include <unordered_set>
#include <string>

#include <iostream>
#include <sstream>
#include <fstream>
#include <istream>

#include<algorithm>
#include <cmath>
#include <typeinfo>


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


namespace io
{

std::vector<std::string>
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

}


namespace basic
{

const std::size_t byteCapacity =
        std::numeric_limits< char >::max() - std::numeric_limits< char >::min();

const std::array< char , 4 > agtc = { 'A' , 'G' , 'T' , 'C' };
const std::array< char , 4 > tcag = { 'T' , 'T' , 'A' , 'G' };
const std::array< int , byteCapacity > codeAGTC([]{
    std::array< int , byteCapacity > codes;
    std::generate_n( codes.begin() , byteCapacity , [](){ return -1;});
    codes['A'] = 0;
    codes['G'] = 1;
    codes['T'] = 2;
    codes['C'] = 3;
    return codes;
}());



/**
 * @brief findClumps
 * Given integers L and t, a string Pattern forms an (L, t)-clump inside a
 * (larger) string Genome if there is an interval of Genome of length L in which
 * Pattern appears at least t times. For example, TGCA forms a (25,3)-clump in the
 * following Genome: gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac.
 *
 * @param input
 * std::string of the sequence to examine the enrichment of kmers in.
 *
 * @param k
 * kmer size of interest.
 *
 * @param windowSize
 * The window size to find the clumps (enrichments) in.
 *
 * @param threshold
 * Occurrances threshold to consider abundant.
 *
 * @return
 * kmers forming clumps as string of vectors.
 */
std::list<std::string>
findClumps( const std::string &input ,
            int k ,
            int windowSize ,
            int threshold )
{
    auto decode = [k]( uint64_t code ){
        std::string s( k , 0 );
        int i=0;
        for( auto i = 0 ; i < k ; i++ ){
            s[ k - i - 1 ] = agtc[ code % 4 ];
            code /= 4;
        }
        return s;
    };

    std::deque< uint64_t > window ;
    std::unordered_set< uint64_t > occurance;
    std::vector< int > occuranceSpace;
    try{
        occuranceSpace = std::vector< int >( std::pow( 4 , k ) , 0 );
    } catch( const std::bad_alloc &e )
    {
        std::cout << e.what();
        exit( EXIT_FAILURE );
    }

    uint64_t sequenceCode = 0;
    const uint64_t mask = static_cast< uint64_t >( std::pow( 4 , k-1 ));
    auto wholeSequence = input.c_str();
    for( auto i = 0; i < k - 1 ; i++ )
        sequenceCode = sequenceCode*4 + codeAGTC[ wholeSequence[ i ]];

    for( auto i = k - 1 ; i < input.size() ; i++ ){
        sequenceCode = ( sequenceCode % mask ) * 4 + codeAGTC[wholeSequence[i]];
        window.push_back( sequenceCode );

        if( i >= windowSize )
        {
            occuranceSpace[ window.front() ]--;
            window.pop_front();
        }

        if( ++occuranceSpace[ sequenceCode ] >= threshold )
        {
            occuranceSpace[ sequenceCode ] = 0;
            occurance.insert( sequenceCode );
        }
    }
    std::list< std::string > frequentWords;
    std::transform( std::begin( occurance ) , std::end( occurance ) ,
                    std::inserter( frequentWords , std::end( frequentWords )) ,
                    decode );
    return frequentWords;
}


std::string
complementSequence( const std::string &sequence )
{
    auto i = sequence.size();
    std::string complement;
    for( auto c : sequence )
        complement[ --i ] = tcag[ codeAGTC[ c ]];
    return complement;
}

std::vector<std::string>
complementSequences( const std::vector<std::string> &sequences )
{
    std::vector< std::string > complemented( sequences.size());
    std::transform( std::begin( sequences ) , std::end( sequences ) ,
                    std::inserter( complemented , std::begin( complemented )) ,
                    complementSequence );
    return complemented;
}


/**
 * @brief skewDiagram
 * Define the skew of a DNA string Genome, denoted Skew(Genome),
 * as the difference between the total number of occurrences of
 * 'G' and 'C' in Genome. Let Prefixi (Genome) denote the prefix
 * (i.e., initial substring) of Genome of length i.
 * @param sequence
 * The genome/sequence to apply the process.
 * @param c1
 * The first term of the gradient = count(c1) - count(c2).
 * @param c2
 * The second term of the gradient = count(c1) - count(c2).
 */
auto
skewDiagram( const std::string &sequence , char c1 = 'G', char c2 = 'C' )
{
    using peak = std::pair< std::list< int > , int >;
    peak min( {} , std::numeric_limits< int >::max());
    peak max( {} , std::numeric_limits< int >::min());
    int skewness = 0;
    for( auto i = 0 ; i < sequence.size() ; i++ )
    {
        if( sequence[ i ] == c1 && ++skewness > max.second )
        {
            max.second = skewness ;
            max.first.clear();
        }
        else if( sequence[ i ] == c2 && --skewness < min.second )
        {
            min.second = skewness ;
            min.first.clear();
        }

        if ( skewness == max.second )
            max.first.push_back( i+1 );

        else if( skewness == min.second )
            min.first.push_back( i+1 );
    }
    return std::make_pair( min , max );
}

/**
 * @brief minimumCoinsChange
 * The Change Problem
 * Find the minimum number of coins needed to make change.
 * @param value
 * An integer money
 * @param domination
 * an array Coins of positive integers.
 * @return The minimum number of coins with denominations Coins that changes money.
 */
int
minimumCoinsChange( int value , const std::list< int > &domination )
{
    std::vector< int > sortedDomination( domination.begin() , domination.end());
    std::sort( sortedDomination.begin() , sortedDomination.end());
    std::vector< int > minCount( value + 1 , std::numeric_limits< int >::max());
    minCount[ 0 ] = 0;
    for( auto money = 1 ; money <= value ; money++ )
    {
        auto count = std::numeric_limits< int >::max();
        for( auto coin : sortedDomination )
        {
            auto remainder = money - coin;
            if( remainder < 0  ) break;
            else if( 1 + minCount[ remainder ] < count )
                    count = 1 + minCount[ remainder ];
        }
        minCount[ money ] = count;
    }
    return minCount[ value ];
}

std::list< std::string >
frequentWordsBruteForce( const std::string &sequence , int kmer , int distance );

int
distance( const std::string &s1 , const std::string &s2 );


}

}
#endif // COMMON_HH
