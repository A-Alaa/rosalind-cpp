#ifndef BA3_HPP
#define BA3_HPP

#include "ba2.hpp"

namespace rosalind
{
namespace ba3
{

/**
 * @brief reconstructStringFromKmers
 * @param firstIt
 * @param lastIt
 * @return
 */
template< typename SeqIt >
std::string
reconstructStringFromKmers( SeqIt firstIt , SeqIt lastIt )
{
    assert( std::distance( firstIt , lastIt ) > 1 );
    std::string str;
    str.reserve( firstIt->length() + std::distance( firstIt + 1 , lastIt ));
    str.insert( str.begin() , firstIt->begin() , firstIt->end());
    std::transform( firstIt + 1 , lastIt ,
                    std::inserter( str , std::end( str )) ,
                    []( const std::string &s ){ return s.back();});
    return str;
}

/**
 * @brief constructOverlapMap
 * @param firstIt
 * @param lastIt
 */
template< typename SeqIt >
auto
constructOverlapMap( SeqIt firstIt , SeqIt lastIt )
{
    using S = std::string;
    using V = std::vector< char >;
    std::unordered_map< S , std::pair< V , V >> prefixSuffixMap;
    // Populate suffixes, O(n).
    for( auto it = firstIt ; it != lastIt ; ++it )
        prefixSuffixMap[ S( it->begin()+1 , it->end() ) ].first.push_back(
                    it->front());

    // Populate prefixes, O(n).
    for( auto it = firstIt ; it != lastIt ; ++it )
        prefixSuffixMap[ S( it->begin() , it->end() - 1 ) ].second.push_back(
                    it->back());

    return prefixSuffixMap;
}

/**
 * @brief overlapStringsGraph
 * Given an arbitrary collection of k-mers Patterns,
 * we form a graph having a node for each k-mer in Patterns
 * and connect k-mers Pattern and Pattern' by a directed edge
 * if Suffix(Pattern) is equal to Prefix(Pattern').
 * The resulting graph is called the overlap graph on
 * these k-mers, denoted Overlap(Patterns).
 * @param firstIt
 * @param lastIt
 * @return
 */
template< typename SeqIt >
std::vector< std::pair< std::string , std::string >>
overlapStringsGraph( SeqIt firstIt , SeqIt lastIt )
{
    // Construct the overlap map, O(n).
    auto prefixSuffixMap = constructOverlapMap( firstIt , lastIt );

    // Construct the graph, O(n).
    std::vector< std::pair< std::string , std::string >> overlap;
    for( const auto &prefixSuffix : prefixSuffixMap )
        for( const auto prefix : prefixSuffix.second.first )
            for( const auto suffix : prefixSuffix.second.second )
            {
                auto prefixStr = prefixSuffix.first;
                auto suffixStr = prefixSuffix.first;
                suffixStr.push_back( suffix );
                overlap.emplace_back(
                            std::make_pair( prefixStr.insert( 0 , 1 , prefix ) ,
                                            suffixStr ));
            }
    return overlap;
}

/**
 * @brief constructDeBruijnGraph
 * @param firstIt
 * @param lastIt
 * @return
 */
template< typename SeqIt >
std::map< std::string , std::vector< std::string >>
constructDeBruijnGraph( SeqIt firstIt , SeqIt lastIt )
{
    // Construct the overlap map, O(n).
    auto prefixSuffixMap = constructOverlapMap( firstIt , lastIt );

    // Construct the graph, O(n).
    std::map< std::string , std::vector< std::string >> deBruijnGraph;
    for( const auto &prefixSuffix : prefixSuffixMap )
        for( const auto prefix : prefixSuffix.second.first )
        {
            auto prefixStr = prefixSuffix.first;
            prefixStr.insert( 0 , 1 , prefix );
            prefixStr.pop_back();
            deBruijnGraph[ prefixStr ].push_back( prefixSuffix.first );
        }
    return deBruijnGraph;
}

/**
 * @brief constructDeBruijnGraph
 * @param input
 */
auto
constructDeBruijnGraph( const RosalindIOType &input )
{
    unsigned int k = atoi( input[0].c_str());
    auto kmers = ba1::extractKmers( input[1] , k );
    return constructDeBruijnGraph( kmers.cbegin() , kmers.cend());
}
}
}
#endif // BA2_HPP
