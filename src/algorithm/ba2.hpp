#ifndef BA2_HPP
#define BA2_HPP

#include "ba1.hpp"

namespace rosalind
{
namespace ba2
{

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

std::vector< std::array< float , 4 >>
makeProfile( const std::vector< std::string > &kmers )
{
    assert( !kmers.empty());
    const unsigned int k = kmers[0].size();
    const unsigned int t = kmers.size();

    assert( std::all_of( kmers.begin() , kmers.end() ,
                         [k]( const std::string &s){ return s.length() == k; }));

    std::vector< std::array< float , 4 >> profile( k , {} );
    for( const auto &kmer : kmers )
        for( IndexType i = 0 ; i < k ; i++ )
            profile[i][ codeACGT[ kmer[ i ]]]++;

    for( auto &histo : profile )
        for( float &value : histo )
            value /= t;

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
std::vector< std::string >
motifEnumeration( const std::vector< std::string > &sequences , unsigned int k , unsigned int d )
{
    std::set< std::string > kmers;
    for( const auto &sequence : sequences )
        for( const auto &kmer : ba1::extractKmers( sequence , k ))
            for( const auto &mutant : ba1::dNeighborhood( kmer , d ))
                kmers.insert( mutant  );

    auto predicate = [k,d]( const std::string &sequence , const std::string &kmer )
    {
        return rosalind::io::findWithMismatches( sequence , kmer , d ) !=
                std::string::npos;
    };

    std::vector< std::string > motifs;
    for( const auto &kmer : kmers )
        if( std::all_of( sequences.cbegin() , sequences.cend() ,
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
    return motifEnumeration( std::vector< std::string >( ++input.cbegin() , input.cend()) ,
                             k , d );
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
    std::vector< std::vector< std::string >> kmers;
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
        CountType minLocalDistance = 0;
        auto distCmp = [&uKmer]( const std::string &a , const std::string &b ){
            return ba1::hammingDistance( uKmer , a ) <
                    ba1::hammingDistance( uKmer , b );
        };
        for( const auto &line : kmers )
            minLocalDistance +=
                    ba1::hammingDistance( uKmer ,
                                          *std::min_element( line.begin() ,
                                                             line.end() ,
                                                             distCmp ));
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
    return medianKmer( std::vector< std::string >( ++input.cbegin() , input.cend()) ,
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
std::vector< std::string >
greedyMotifSearch( const std::vector< std::string > &sequences ,
                   unsigned int k )
{
    const unsigned int t = sequences.size();
    assert( t > 1 );
    std::vector< std::vector< std::string >> kmers;
    std::transform( std::begin( sequences ) , std::end( sequences ) ,
                    std::inserter( kmers , std::end( kmers )) ,
                    [k]( const std::string &sequence ){ return ba1::extractKmers( sequence , k );});

    auto score = []( const std::vector< std::string > &kmers )
    {
        auto consensusMotif = consensus( kmers );
        return std::accumulate( std::begin( kmers ) , std::end( kmers ) ,
                                0 , [consensusMotif]( unsigned int a , const std::string &kmer ){
            return a + ba1::hammingDistance( consensusMotif , kmer );
        });
    };

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
            auto profile = makeProfile( V( motifs.begin() ,
                                           motifs.begin() + j ));
            motifs.push_back( profileMostProbableKmer( sequences[ j ] , profile ));
        }
        unsigned int motifsScore = score( motifs );
        if( motifsScore < bestScore )
        {
            bestScore = motifsScore;
            bestMotifs = motifs;
        }
    }
    return bestMotifs;
}

/**
 * @brief greedyMotifSearch
 * @param input
 */
auto
greedyMotifSearch( const RosalindIOType &input )
{
    const unsigned int k = std::atoi( io::split( input[0] , ' ')[0].c_str());
    return greedyMotifSearch( std::vector< std::string >( input.begin() + 1 ,
                                                          input.end()) , k );
}
}
}
#endif // BA2_HPP
