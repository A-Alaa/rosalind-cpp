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


namespace basic
{

const std::size_t byteCapacity =
        std::numeric_limits< char >::max() - std::numeric_limits< char >::min();

const std::array< char , 4 > acgt = { 'A' , 'C' , 'G' , 'T' };
const std::array< char , 4 > tgca = { 'T' , 'G' , 'C' , 'A' };
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

/**
 * @brief extractKmers
 * @param text
 * @param k
 * @return
 */
std::vector< std::string > extractKmers( const std::string &text , int k )
{
    std::vector< std::string > kmers;
    int n = text.length();
    for ( auto i = 0 ; i < n-k+1 ; i++ )
        kmers.push_back( text.substr( i , k ));
    return kmers;
}

/**
 * ba1m
 * @brief decode
 * Convert an integer to its corresponding DNA string.
 * @param code
 * @param k
 * @return
 * NumberToPattern(index, k).
 */
template< typename T >
std::string decode( T code , int k ){
    std::string s( k , 0 );
    for( auto i = 0 ; i < k ; i++ )
    {
        s[ k - i - 1 ] = acgt[ code % 4 ];
        code /= 4;
    }
    return s;
}

/**
 * @brief decode
 * @param input
 */
template< typename T = CodeType >
auto numberToPattern( const RosalindIOType &input ){
    return decode< T >(
                std::atoll( input[0].c_str( )) , std::atoi( input[1].c_str( )));
}

/**
 * ba1l
 * @brief encode
 * Convert a DNA string to a number.
 * @param pattern
 * @return
 * PatternToNumber(Pattern).
 */
template< typename T = CodeType >
T encode( const std::string &pattern )
{
    T code = 0;
    for( const auto c : pattern )
        code = code * 4 + codeACGT[ c ];
    return code;
}

/**
 * @brief encode
 * @param input
 */
auto
encode( const RosalindIOType &input )
{
    return encode< uint64_t >( input[0] );
}


/**
 * ba1a
 * @brief patternCount
 * @param sequence
 * @param pattern
 * @return
 */
auto
patternCount( const std::string &sequence , const std::string &pattern )
{
    CountType count = 0 ;
    IndexType nextIdx = sequence.find( pattern );
    while( nextIdx != std::string::npos )
    {
        count++;
        nextIdx = sequence.find( pattern , nextIdx + 1 );
    }
    return count;
}

auto
patternCount( const RosalindIOType &input )
{
    return patternCount( input[0] , input[1] );
}

/**
 * ba1b
 * @brief frequentWordsBruteForce
 * Find the most frequent k-mers in a string.
 * We say that Pattern is a most frequent k-mer in Text
 * if it maximizes Count(Text, Pattern) among all k-mers.
 * For example, "ACTAT" is a most frequent 5-mer in "ACAACTATGCATCACTATCGGGAACTATCCT",
 * and "ATA" is a most frequent 3-mer of "CGATATATCCATAG".
 * @param sequence
 * A string to look for frequent words in.
 * @param k
 * Kmers fixed size of interest.
 * @return
 * Most frequent kmers.
 */
std::list< std::string >
frequentWordsBruteForce( const std::string &sequence , int k )
{
    std::unordered_map< std::string , int > frequency;
    const auto kmers = extractKmers( sequence , k );
    for( auto &kmer : kmers )
        frequency[ kmer ]++;

    std::list< std::string > mostFrequentKmers;
    int maxFrequency = 0;
    for( auto &kmer : frequency )
    {
        if( kmer.second > maxFrequency )
        {
            mostFrequentKmers.clear();
            maxFrequency = kmer.second;
        }
        if( kmer.second == maxFrequency )
            mostFrequentKmers.push_back( kmer.first );
    }
    return mostFrequentKmers;
}

/**
 * @brief frequentWordsBruteForce
 * @param sequence
 */
auto
frequentWordsBruteForce( const RosalindIOType &inputStrings )
{
    return frequentWordsBruteForce( inputStrings[ 0 ] ,
            std::atoi( inputStrings[ 1 ].c_str( )));
}

/**
 * ba1c
 * @brief complementSequence
 * Find the reverse complement of a DNA string.
 *
 * In DNA strings, symbols 'A' and 'T' are complements of each other,
 * as are 'C' and 'G'. Given a nucleotide p, we denote its complementary
 * nucleotide as p. The reverse complement of a DNA string
 * Pattern = p1…pn is the string Pattern = pn … p1 formed by
 * taking the complement of each nucleotide in Pattern,
 * then reversing the resulting string.
 * For example, the reverse complement of
 * Pattern = "GTCA" is Pattern = "TGAC".
 * @param sequence
 * The sequence to find its complement.
 * @return
 * The complemented sequence.
 */
std::string
complementSequence( const std::string &sequence )
{
    auto i = sequence.size();
    std::string complement( sequence.size() , 0 );
    for( auto c : sequence )
        complement[ --i ] = tgca[ codeACGT[ c ]];
    return complement;
}


/**
 * ba1c
 * @brief complementSequence2
 * @param sequence
 * @return
 */
std::string
complementSequence2( const std::string &sequence )
{
    std::string complement;
    std::transform( sequence.rbegin() , sequence.rend() ,
                    std::back_inserter( complement ) ,
                    []( char c ){ return tgca[ codeACGT[ c ]];});
    return complement;
}


std::vector<std::string>
complementSequences( const std::vector<std::string> &sequences )
{
    std::vector< std::string > complemented;
    std::transform( std::begin( sequences ) , std::end( sequences ) ,
                    std::inserter( complemented , std::begin( complemented )) ,
                    complementSequence );
    return complemented;
}

/**
 * ba1d
 * @brief patternMatching
 * Find all occurrences of a pattern in a string.
 * @param sequence
 * The big sequence where to look for the pattern.
 * @param pattern
 * The pattern to search in sequence.
 * @return
 * vector of indices corresponds to occurances of pattern in sequence.
 */
std::vector< IndexType >
patternMatching( const std::string &sequence , const std::string &pattern )
{
    std::vector< IndexType > occurances;
    auto nextIndex = sequence.find( pattern );
    while( nextIndex != std::string::npos )
    {
        occurances.push_back( nextIndex );
        nextIndex = sequence.find( pattern , nextIndex + 1 );
    }
    return occurances;
}

/**
 * @brief patternMatching
 * @param input
 */
auto
patternMatching( const RosalindIOType &input )
{
    return patternMatching( input[1] , input[0] );
}

/**
 * BA1E
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
std::list< std::string >
findClumps( const std::string &input ,
            unsigned int k ,
            unsigned int windowSize ,
            unsigned int threshold )
{
    std::deque< IndexType > window ;
    std::unordered_set< IndexType > occurance;
    std::vector< CountType > occuranceSpace;
    try{
        occuranceSpace = std::vector< CountType >( powi( 4 , k ) , 0 );
    } catch( const std::bad_alloc &e )
    {
        std::cout << e.what();
        exit( EXIT_FAILURE );
    }

    CodeType sequenceCode = 0;
    const CodeType mask = powi( 4 , k-1 );
    auto wholeSequence = input.c_str();

    for( IndexType i = 0; i < k - 1 ; i++ )
        sequenceCode = sequenceCode*4 + codeACGT[ wholeSequence[ i ]];

    for( IndexType i = k - 1 ; i < input.size() ; i++ )
    {
        sequenceCode = ( sequenceCode % mask ) * 4 + codeACGT[wholeSequence[i]];
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
                    std::bind( decode< CodeType > , std::placeholders::_1 , k ));
    return frequentWords;
}

/**
 * @brief findClumps
 * @param input
 */
auto
findClumps( const RosalindIOType &input )
{
    auto parameters = rosalind::io::split( input[1] , ' ');
    assert( parameters.size() == 3 );
    return findClumps( input[ 0 ] ,
            atoi( parameters[0].c_str( )) ,
            atoi( parameters[1].c_str( )) ,
            atoi( parameters[2].c_str( )));
}

/**
 * BA1F
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
    using peak = std::pair< std::list< IndexType > , int >;
    peak min( {} , std::numeric_limits< int >::max());
    peak max( {} , std::numeric_limits< int >::min());
    int skewness = 0;
    for( IndexType i = 0 ; i < sequence.size() ; i++ )
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
 * ba1g
 * @brief hammingDistance
 * Compute the Hamming distance between two DNA strings.
 * We say that position i in k-mers p1 … pk and q1 … qk is
 * a mismatch if pi ≠ qi. For example, CGAAT and CGGAC have
 *  two mismatches. The number of mismatches between
 * strings p and q is called the Hamming distance between
 * these strings and is denoted HammingDistance(p, q).
 * @param s1
 * @param s2
 * @return
 * An integer value representing the Hamming distance.
 */
CountType
hammingDistance( const std::string &s1 , const std::string &s2 )
{
    const auto size = ( s1.size() < s2.size())? s1.size() : s2.size();
    CountType distance = 0;
    for( IndexType i = 0 ; i < size ; i++ )
        distance += s1[ i ] != s2[ i ];
    return distance;
}

/**
 * ba1g
 * @brief hammingDistance
 * @param s1
 * @param s2
 * @return
 */
CountType
hammingDistance( const char * const s1 , const char * const s2 , size_t size )
{
    CountType distance = 0;
    for( IndexType i = 0 ; i < size ; i++ )
        distance += s1[ i ] != s2[ i ];
    return distance;
}


/**
 * ba1h
 * @brief approximatePatternMatching
 * Find all approximate occurrences of a pattern in a string.
 * We say that a k-mer Pattern appears as a substring of Text
 * with at most d mismatches if there is some k-mer substring
 * Pattern' of Text having d or fewer mismatches with Pattern,
 * i.e., HammingDistance(Pattern, Pattern') ≤ d. Our observation
 * that a DnaA box may appear with slight variations leads to the
 * following generalization of the Pattern Matching Problem.
 * @param sequence
 * @param pattern
 * @param distance
 * @return
 * All starting positions where Pattern appears as a substring of Text with at most d mismatches.
 */
std::vector< IndexType >
approximatePatternMatching( const std::string &sequence ,
                            const std::string &pattern ,
                            CountType distance = 0 )
{
    assert( sequence.size() >= pattern.size());
    std::vector< IndexType > indices;
    const auto size = pattern.size();
    const IndexType indexSpace = sequence.size() - pattern.size();
    auto cSequence = sequence.c_str();
    auto cPattern  = pattern .c_str();
    for( IndexType i = 0 ; i < indexSpace + 1 ; i++ )
        if( hammingDistance( cSequence + i , cPattern , size ) <= distance )
            indices.push_back( i );
    return indices;
}

auto
approximatePatternMatching( const RosalindIOType &input )
{
    return rosalind::basic::approximatePatternMatching(
                input[1] , input[0] ,  atoi( input[2].c_str() ));
}

/**
 * ba1i
 * @brief frequentWordsWithMismatches
 * Find the most frequent k-mers with mismatches in a string.
 * Given strings Text and Pattern as well as an integer d,
 * we define Countd(Text, Pattern) as the total number of
 * occurrences of Pattern in Text with at most d mismatches.
 * For example, Count1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4
 * because AAAAA appears four times in this string with at most
 * one mismatch: AACAA, ATAAA, AAACA, and AAAGA. Note that two of
 * these occurrences overlap.
 * A most frequent k-mer with up to d mismatches in Text
 * is simply a string Pattern maximizing Countd(Text, Pattern)
 * among all k-mers. Note that Pattern does not need to actually
 * appear as a substring of Text; for example, AAAAA is the most
 * frequent 5-mer with 1 mismatch in AACAAGCTGATAAACATTTAAAGAG,
 * even though AAAAA does not appear exactly in this string.
 *
 * @param sequence
 * @param k
 * Kmers fixed size.
 * @param d
 * Maximum distance to consider an approximate match.
 * @return
 * All most frequent k-mers with up to d mismatches in sequence.
 */
std::list< std::string >
frequentWordsWithMismatches( const std::string &sequence ,
                             unsigned int k , unsigned int d )
{
    auto cSequence = sequence.c_str();
    std::pair< std::list< std::string > , int > mostFrequentKmers;
    const CodeType kmerSpace = static_cast< CodeType >( std::pow( 4 , k ));
    for( IndexType i = 0 ; i < kmerSpace ; i++ )
    {
        int occurance = 0;
        auto kmer = decode( i , k );
        auto cKmer = kmer.c_str();
        for( IndexType j = 0 ; j < sequence.size() - k + 1 ; j++ )
            occurance += hammingDistance( cSequence + j , cKmer , k ) <= d;

        if( occurance > mostFrequentKmers.second )
        {
            mostFrequentKmers.second = occurance;
            mostFrequentKmers.first.clear();
        }
        if( occurance == mostFrequentKmers.second )
            mostFrequentKmers.first.push_back( kmer );
    }
    return mostFrequentKmers.first;
}

auto
frequentWordsWithMismatches( const RosalindIOType &inputStrings  )
{
    auto parameters = rosalind::io::split( inputStrings[1] , ' ');
    int k = atoi( parameters[0].c_str());
    int d = atoi( parameters[1].c_str());
    return frequentWordsWithMismatches( inputStrings[0] , k , d );
}


/**
 * ba1j
 * @brief frequentWordsWithMismatchesAndReverseComplement
 * Find the most frequent k-mers (with mismatches and reverse complements) in a DNA string.
 * @param sequence
 * @param k
 * @param d
 * @return
 * All k-mers Pattern maximizing the
 * sum Countd(Text, Pattern) + Countd(Text, Pattern)
 * over all possible k-mers.
 */
std::list< std::string >
frequentWordsWithMismatchesAndReverseComplement(
        const std::string &sequence ,
        unsigned int k , unsigned int d )
{
    auto cSequence = sequence.c_str();
    std::pair< std::list< std::string > , int > mostFrequentKmers;
    const CodeType kmerSpace = powi( 4 , k );
    for( IndexType i = 0 ; i < kmerSpace ; i++ )
    {
        int occurance = 0;
        auto kmer = decode( i , k );
        auto rKmer = rosalind::basic::complementSequence( kmer );
        auto cKmer = kmer.c_str();
        auto crKmer = rKmer.c_str();

        for( IndexType j = 0 ; j < sequence.size() - k + 1 ; j++ )
        {
            occurance += hammingDistance( cSequence + j , cKmer , k ) <= d;
            occurance += hammingDistance( cSequence + j , crKmer , k ) <= d ;
        }
        if( occurance > mostFrequentKmers.second )
        {
            mostFrequentKmers.second = occurance;
            mostFrequentKmers.first.clear();
        }
        if( occurance == mostFrequentKmers.second )
            mostFrequentKmers.first.push_back( kmer );
    }
    return mostFrequentKmers.first;
}

/**
 * @brief frequentWordsWithMismatches
 * @param inputStrings
 */
auto
frequentWordsWithMismatchesAndReverseComplement( const RosalindIOType &inputStrings  )
{
    auto parameters = rosalind::io::split( inputStrings[1] , ' ');
    int k = atoi( parameters[0].c_str());
    int d = atoi( parameters[1].c_str());
    return frequentWordsWithMismatchesAndReverseComplement( inputStrings[0] , k , d );
}

/**
 * @brief stringFrequencyArray
 * Given an integer k, we define the frequency array of a string
 * Text as an array of length 4k, where the i-th element of the array
 * holds the number of times that the i-th k-mer (in the lexicographic order)
 * appears in Text
 * @param sequence
 * @param k
 * @return
 * The frequency array of k-mers in Text.
 */
std::vector< CountType >
stringFrequencyArray( const std::string & sequence , unsigned int k )
{
    auto kmerSpace = powi( 4 , k );
    std::vector< CountType > frequencyArray( kmerSpace , 0 );
    for( const auto &kmer : extractKmers( sequence , k ))
        frequencyArray[ encode( kmer )]++;
    return frequencyArray;
}

/**
 * @brief stringFrequencyArray
 * @param input
 */
auto
stringFrequencyArray( const RosalindIOType &input )
{
    return stringFrequencyArray( input[0] , atoi( input[1].c_str( )));
}

/**
 * BA5A
 * @brief minimumCoinsChange
 * The Change Problem
 * Find the minimum number of coins needed to make change.
 * @param value
 * An integer money
 * @param domination
 * an array Coins of sorted positive integers.
 * @return The minimum number of coins with denominations Coins that changes money.
 */
auto
minimumCoinsChange( int value , const std::vector< int > &domination )
{
    std::vector< CountType > minCount( value + 1 ,
                                       std::numeric_limits< CountType >::max());
    minCount[ 0 ] = 0;
    for( auto money = 1 ; money <= value ; money++ )
    {
        auto min = std::numeric_limits< CountType >::max();
        for( auto coin : domination )
        {
            auto remainder = money - coin;
            if( remainder < 0  ) break;
            else if( 1 + minCount[ remainder ] < min )
                min = 1 + minCount[ remainder ];
        }
        minCount[ money ] = min;
    }
    return minCount[ value ];
}

/**
 * BA5A
 * @brief minimumCoinsChange
 * @param money
 * @param domination
 * @param delim
 * @return
 */
auto
minimumCoinsChange( const RosalindIOType &input )
{
    auto money = std::atoi( input[0].c_str() );
    auto _domination = rosalind::io::split( input[1] , ',' );
    std::vector< int > __domination;
    std::transform( std::begin( _domination ) , std::end( _domination ) ,
                    std::inserter( __domination , std::begin( __domination )) ,
                    []( const std::string &s ){ return std::atoi( s.c_str());} );
    std::sort( __domination.begin() , __domination.end());
    return minimumCoinsChange( money , __domination );
}




}

}
#endif // COMMON_HH
