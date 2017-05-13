#ifndef BA2_HPP
#define BA2_HPP

#include "ba1.hpp"

namespace rosalind
{
namespace ba2
{

/**
 * @brief mapKmers
 * @param sequences
 * @return
 */
template< typename SeqIt >
std::vector< std::vector< std::string >>
mapKmers( SeqIt firstIt , SeqIt lastIt , unsigned int k )
{
    std::vector< std::vector< std::string >> kmers;
    std::transform( firstIt , lastIt ,
                    std::inserter( kmers , std::end( kmers )) ,
                    [k]( const std::string &sequence ){ return ba1::extractKmers( sequence , k );});
    return kmers;
}

/**
 * @brief consensus
 * @param kmers
 * @return
 */
std::string
consensus( const std::vector< std::string > &kmers )
{
    assert( !kmers.empty());
    const unsigned int k = kmers[0].size();
    const unsigned int t = kmers.size();
    assert( std::all_of( kmers.begin() , kmers.end() ,
                         [k]( const std::string &s){ return s.length() == k; }));

    std::vector< std::array< CountType , 4 >> frequency( k , {});
    for( const auto &kmer : kmers )
        for( unsigned int i = 0 ; i < k ; i++ )
            frequency[i][ codeACGT[ kmer[ i ]]]++;
    std::string consensus;
    std::transform( std::begin( frequency ) , std::end( frequency ) ,
                    std::inserter( consensus , consensus.end()) ,
                    []( const std::array< CountType , 4 > &histo )
    {
        return acgt[ std::max_element( histo.begin() , histo.end()) - histo.begin()];
    });
    return consensus;
}

/**
 * @brief consensus
 * @param profile
 * @return
 */
std::string
consensus( const std::vector< std::array< CountType , 4 >> &profile )
{
    assert( !profile.empty());
    std::string consensus;
    std::transform( std::begin( profile ) , std::end( profile ) ,
                    std::inserter( consensus , consensus.end()) ,
                    []( const std::array< CountType , 4 > &histo )
    {
        return acgt[ std::max_element( histo.begin() , histo.end()) - histo.begin()];
    });
    return consensus;
}

/**
 * @brief hammingDistanceScore
 * @param kmers
 * @return
 */
int
hammingDistanceScore( const std::vector< std::string > &kmers )
{
    auto consensusMotif = consensus( kmers );
    return std::accumulate( std::begin( kmers ) , std::end( kmers ) ,
                            0 , [consensusMotif]( int a , const std::string &kmer ){
        return a + ba1::hammingDistance( consensusMotif , kmer );
    });
}

/**
 * @brief makeProfile
 * @param kmers
 * @param pseudoCount
 * @return
 */
template< typename SeqIt >
std::vector< std::array< float , 4 >>
makeProfile( SeqIt firstIt , SeqIt lastIt ,
             float pseudoCount = 0 )
{
    const unsigned int k = (*firstIt).length();
    const unsigned int t = std::distance( firstIt , lastIt);
    assert( t > 0 );


    assert( std::all_of( firstIt , lastIt ,
                         [k]( const std::string &s){ return s.length() == k; }));

    std::vector< std::array< float , 4 >> profile( k , {} );
    for( auto it = firstIt ; it != lastIt ; it++ )
        for( IndexType i = 0 ; i < k ; ++i )
            profile[i][ codeACGT[ (*it)[ i ]]]++;

    for( auto &histo : profile )
        for( float &value : histo )
            value = ( value + pseudoCount ) / t;

    return profile;
}

template< typename SeqIt >
std::vector< std::array< float , 4 >>
makeProfileWithException( SeqIt firstIt , SeqIt lastIt ,
                          std::vector< std::string >::const_iterator exceptIt ,
                          float pseudoCount = 0 )
{
    const unsigned int k = (*firstIt).length();
    const unsigned int t = std::distance( firstIt , lastIt );
    assert( t > 1 );

    assert( std::all_of( firstIt , lastIt ,
                         [k]( const std::string &s){ return s.length() == k; }));

    std::vector< std::array< float , 4 >> profile( k , {} );
    for( auto it = firstIt ; it != lastIt ; ++it )
        if ( it == exceptIt ) continue;
        else
            for( IndexType i = 0 ; i < k ; ++i )
                profile[i][ codeACGT[ (*it)[ i ]]]++;

    for( auto &histo : profile )
        for( float &value : histo )
            value = ( value + pseudoCount ) / ( t - 1 );

    return profile;
}

/**
 * ba2a
 * @brief motifEnumeration
 * Given a collection of strings Dna and an integer d,
 * a k-mer is a (k,d)-motif if it appears in every
 * string from Dna with at most d mismatches.
 * @param sequences
 * @param k
 * @param d
 * @return
 * All (k, d)-motifs in Dna.
 */
template< typename SeqIt >
std::vector< std::string >
motifEnumeration( SeqIt firstIt , SeqIt lastIt ,
                  unsigned int k , unsigned int d )
{
    std::set< std::string > kmers;
    for( auto seqIt = firstIt ; seqIt != lastIt ; ++seqIt )
        for( const auto &kmer : ba1::extractKmers( *seqIt , k ))
            for( const auto &mutant : ba1::dNeighborhood( kmer , d ))
                kmers.insert( mutant  );

    auto predicate = [k,d]( const std::string &sequence , const std::string &kmer )
    {
        return rosalind::io::findWithMismatches( sequence , kmer , d ) !=
                std::string::npos;
    };

    std::vector< std::string > motifs;
    for( const auto &kmer : kmers )
        if( std::all_of( firstIt , lastIt ,
                         std::bind( predicate , std::placeholders::_1 , kmer )))
            motifs.push_back( kmer );

    return motifs;
}

/**
 * @brief motifEnumeration
 * @param input
 */
auto
motifEnumeration( const RosalindIOType &input )
{
    auto parameters = rosalind::io::split( input[0] , ' ');
    auto k = std::atoi( parameters[0].c_str());
    auto d = std::atoi( parameters[1].c_str());
    return motifEnumeration( input.cbegin() + 1 , input.cend() ,
                             k , d );
}


/**
 * ba2h
 * @brief minimalHammingDistance
 * The first potential issue with implementing MedianString
 * from “Find a Median String” is writing a function to compute
 * d(Pattern, Dna) = the sum of distances between Pattern
 * and each string in Dna = {Dna1, ..., Dnat}.
 * @param pattern
 * pattern of interest to measure the minimal hamming distance accross
 * each string.
 * @param firstIt
 * Iterator pointing to the first string sequence in container.
 * @param lastIt
 * Iterator pointing to the last (excluded) string in a container.
 * @return
 */
template< typename SeqIt >
CountType
minimalHammingDistance( const std::string &pattern ,
                        SeqIt firstIt , SeqIt lastIt )
{
    //    static_assert( std::is_same< SeqIt , std::vector< std::vector< std::string >>::iterator >::value ,
    //                   "SeqIt must be an iterator to containers of std::string(s)");
    CountType minimalDistance = 0;
    auto distCmp = [&pattern]( const std::string &a , const std::string &b ){
        return ba1::hammingDistance( pattern , a ) <
                ba1::hammingDistance( pattern , b );
    };
    using K = std::vector< std::vector< std::string >>;
    for( K::iterator lineIt = firstIt ; lineIt != lastIt ; ++lineIt )
        minimalDistance +=
                ba1::hammingDistance( pattern ,
                                      *std::min_element( lineIt->begin() ,
                                                         lineIt->end() ,
                                                         distCmp ));
    return minimalDistance;
}

auto
minimalHammingDistance( const RosalindIOType &input )
{
    auto kmer = input[0];
    auto sequences = io::split( input[1] , ' ');
    const unsigned int k = kmer.length();

    using KmersType = std::vector< std::vector< std::string >>;
    KmersType kmers;
    std::transform( std::begin( sequences ) , std::end( sequences ),
                    std::inserter( kmers , kmers.end()) ,
                    [k]( const std::string &sequence ){
        return ba1::extractKmers( sequence , k ) ;});

    return minimalHammingDistance( input[0] , kmers.begin() , kmers.end());
}

/**
 * ba2b
 * @brief medianKmer
 * Find a median string.
 * @param sequences
 * @param k
 * @return
 * A k-mer Pattern that minimizes
 * distance(Pattern, Dna) over all k-mers Pattern.
 * (If multiple answers exist, you may return any one.)
 */
std::string medianKmer( const std::vector< std::string > &sequences , unsigned int k )
{
    using KmersType = std::vector< std::vector< std::string >>;
    KmersType kmers;
    std::set< std::string > uniqueKmers;

    std::transform(
                std::begin( sequences ) , std::end( sequences ),
                std::inserter( kmers , kmers.end()) ,
                [&uniqueKmers,k]( const std::string &sequence ){
        auto _ = ba1::extractKmers( sequence , k );
        uniqueKmers.insert( _.begin() , _.end());
        return _ ;});

    std::string medianKmer;
    CountType minDistance = std::numeric_limits< CountType >::max();
    for( const auto &uKmer : uniqueKmers )
    {
        auto minLocalDistance =
                minimalHammingDistance( uKmer , kmers.begin() , kmers.end());
        if( minLocalDistance < minDistance )
        {
            minDistance = minLocalDistance;
            medianKmer = uKmer;
        }
    }
    return medianKmer;
}

/**
 * @brief medianKmer
 * @param input
 */
auto
medianKmer( const RosalindIOType &input )
{
    return medianKmer( std::vector< std::string >( input.cbegin() + 1 , input.cend()) ,
                       std::atoi( input[0].c_str()));
}

/**
 * ba2c
 * @brief profileMostProbableKmer
 * Find a Profile-most probable k-mer in a string.
 * @param sequence
 * @return
 * A Profile-most probable k-mer in Text.
 */
std::string
profileMostProbableKmer( const std::string &sequence ,
                         const std::vector< std::array< float , 4 >> &profile )
{
    const unsigned int L = sequence.size();
    const unsigned int k = profile.size();
    const auto kmers = ba1::extractKmers( sequence , k );

    assert( k <= L );
    auto convolute = [k,&profile]( const std::string &kmer )
    {
        float product = 1.f;
        for( IndexType i = 0 ; i < k ; i++ )
            product *= profile.at( i ).at( codeACGT[ kmer[ i ]] );
        return product;
    };

    IndexType maxIndex = 0;
    float maxProbability = 0;
    for( IndexType i = 0 ; i <= L - k ; i++ )
    {
        float probability = convolute( kmers[i] );
        if( probability > maxProbability )
        {
            maxProbability = probability;
            maxIndex = i;
        }
    }
    return kmers[ maxIndex ];
}

/**
 * @brief profileMostProbableKmer
 * @param input
 */
auto
profileMostProbableKmer( const RosalindIOType &input )
{
    assert( input.size() == 6 );

    auto sequence = input[0];
    auto k = std::atoi( input[1].c_str());
    auto a = io::split( input[2] , ' ');
    auto c = io::split( input[3] , ' ');
    auto g = io::split( input[4] , ' ');
    auto t = io::split( input[5] , ' ');
    std::vector< std::array< float , 4 >> profile;

    for( auto i = 0 ; i < k ; i++ )
        profile.push_back( { std::strtof( a[i].c_str() , 0 ),
                             std::strtof( c[i].c_str() , 0 ),
                             std::strtof( g[i].c_str() , 0 ),
                             std::strtof( t[i].c_str() , 0 )});

    return profileMostProbableKmer( sequence , profile );
}

/**
 * ba2d
 * @brief greedyMotifSearch
 * @param sequences
 * @param k
 * @return
 */
template< typename SeqIt >
std::vector< std::string >
greedyMotifSearch( SeqIt firstIt , SeqIt lastIt ,
                   unsigned int k ,
                   float pseudoCount = 0 )
{
    const unsigned int t = std::distance( firstIt , lastIt );
    assert( t > 1 );
    auto kmers = mapKmers( firstIt , lastIt , k );

    unsigned int bestScore = std::numeric_limits< int >::max();
    std::vector< std::string > bestMotifs;
    for( const auto &line : kmers )
        bestMotifs.push_back( line[0] );

    for( const auto &motif : kmers[0] )
    {
        decltype( bestMotifs ) motifs;
        motifs.push_back( motif );
        for( unsigned int j = 1 ; j < t ; j++ )
        {
            using V = std::vector< std::string >;
            auto profile = makeProfile( motifs.cbegin() ,
                                        motifs.cbegin() + j  ,
                                        pseudoCount );
            motifs.push_back( profileMostProbableKmer( *( firstIt + j ) , profile ));
        }
        unsigned int motifsScore = hammingDistanceScore( motifs );
        if( motifsScore < bestScore )
        {
            bestScore = motifsScore;
            bestMotifs = motifs;
        }
    }
    return bestMotifs;
}

/**
 *
 * @brief greedyMotifSearch
 * @param input
 */
auto
greedyMotifSearch( const RosalindIOType &input )
{
    const unsigned int k = std::atoi( io::split( input[0] , ' ')[0].c_str());
    return greedyMotifSearch( input.cbegin() + 1 , input.cend() , k );
}

/**
 * ba2e
 * @brief greedyMotifSearchWithPseudoCount
 * @param input
 */
auto
greedyMotifSearchWithPseudoCount( const RosalindIOType &input )
{
    const unsigned int k = std::atoi( io::split( input[0] , ' ')[0].c_str());
    return greedyMotifSearch( input.begin() + 1 , input.cend() , k , 1.f );
}

/**
 * ba2f
 * @brief randomizedMotifSearch
 * @param sequences
 * @param k
 * @param pseudoCount
 * @return
 * A collection BestMotifs resulting from running RandomizedMotifSearch(Dna, k, t)
 * 1000 times. Remember to use pseudocounts!
 */
template< typename SeqIt >
std::vector< std::string >
randomizedMotifSearch( SeqIt firstIt , SeqIt lastIt ,
                       unsigned int k ,
                       unsigned int runs = 1000,
                       float pseudoCount = 1.f )
{
    const unsigned int t = std::distance( firstIt , lastIt );
    assert( t > 1 );
    auto kmers = mapKmers( firstIt , lastIt , k );

    auto bestScore = std::numeric_limits< int >::max();
    using MotifsType = std::vector< std::string > ;
    MotifsType bestMotifs;
    while( runs-- > 0 )
    {
        MotifsType oldMotifs;
        std::transform( std::begin( kmers ) , std::end( kmers ) ,
                        std::inserter( oldMotifs , oldMotifs.end()) ,
                        []( const MotifsType &line ){
            return *randomElement( line.begin() , line.end());
        });

        bool converging = true;
        while( converging )
        {
            auto profile = makeProfile( oldMotifs.cbegin() ,
                                        oldMotifs.cend() ,
                                        pseudoCount );
            MotifsType newMotifs;
            std::transform( firstIt , lastIt ,
                            std::inserter( newMotifs , newMotifs.end()) ,
                            [&profile]( const std::string &sequence ){
                return profileMostProbableKmer( sequence , profile );
            });
            auto oldScore = hammingDistanceScore( oldMotifs );
            auto newScore = hammingDistanceScore( newMotifs );
            if(( converging = newScore < oldScore ))
                oldMotifs = newMotifs;
        }

        auto currentScore = hammingDistanceScore( oldMotifs );
        if( currentScore < bestScore )
        {
            bestScore = currentScore;
            bestMotifs = oldMotifs;
        }
    }
    return bestMotifs;
}

/**
 * @brief randomizedMotifSearch
 * @param input
 */
auto
randomizedMotifSearch( const RosalindIOType &input )
{
    const unsigned int k = std::atoi( io::split( input[0] , ' ')[0].c_str());
    return randomizedMotifSearch( input.cbegin() + 1 , input.cend() , k , 1000 , 1.f );
}

/**
 * ba2g
 * @brief gibbsSampler
 * @param firstIt
 * @param lastIt
 * @param k
 * @param runs
 * @param randomStarts
 * @return
 */
template< typename SeqIt >
std::vector< std::string >
gibbsSampler( SeqIt firstIt , SeqIt lastIt , unsigned int k ,
              const int runs , const int randomStarts = 60 )
{
    const unsigned int t = std::distance( firstIt , lastIt );
    assert( t > 1 );
    auto kmers = mapKmers( firstIt , lastIt , k );

    using MotifsType = std::vector< std::string > ;
    MotifsType bestMotifs;
    auto bestScore = std::numeric_limits< int >::max();
    auto _randomStarts = randomStarts;
    while( _randomStarts-- > 0 )
    {
        MotifsType motifs;
        std::transform( std::begin( kmers ) , std::end( kmers ) ,
                        std::inserter( motifs , motifs.end()) ,
                        []( const MotifsType &line ){
            return *randomElement( line.begin() , line.end());
        });

        auto currentBestScore = hammingDistanceScore( motifs );
        MotifsType currentBestMotifs = motifs;
        auto _runs = runs;
        while( _runs-- > 0 )
        {
            auto motifIt = randomElement( motifs.begin() , motifs.end());
            auto idx = std::distance( motifs.begin() , motifIt );
            auto profile = makeProfileWithException( motifs.cbegin() ,
                                                     motifs.cend() , motifIt , 1.f );
            *motifIt = profileMostProbableKmer( *(firstIt + idx ) , profile );
            auto score = hammingDistanceScore( motifs );
            if( score < currentBestScore )
            {
                currentBestMotifs = motifs;
                currentBestScore = score;
            }
        }

        if( currentBestScore < bestScore )
        {
            bestScore = currentBestScore;
            bestMotifs = currentBestMotifs;
        }
    }
    return bestMotifs;
}

/**
 * @brief gibbsSampler
 * @param input
 */
auto
gibbsSampler( const RosalindIOType &input )
{
    auto parameters = io::split( input[0] , ' ');
    unsigned int k = atoi( parameters[0].c_str());
    int n = atoi( parameters[2].c_str());
    return gibbsSampler( input.cbegin() + 1 , input.cend() , k , 2 * n );
}
}
}
#endif // BA2_HPP
