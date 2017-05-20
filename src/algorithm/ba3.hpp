#ifndef BA3_HPP
#define BA3_HPP

#include "ba2.hpp"

namespace rosalind
{
namespace ba3
{

template< typename V >
class Graph
{
public:
    class Vertex
    {
        friend class Graph;
        V _data;
    public:
        Vertex( V data ) : _data( data ){}
        Vertex ( const Vertex & ) = default;
        Vertex() = default;
        const V &data() const
        {
            return _data;
        }
        bool operator ==( const Vertex &other ) const
        {
            _data == other._data;
        }
    };

    class Edge
    {
        friend class Graph;
        Vertex _src;
        Vertex _trgt;
    public:
        Edge( const Vertex &src , const Vertex &target )
            : _src( src ), _trgt( target ){}
        Edge ( const Edge & ) = default;
        Edge() = default ;

        bool operator==( const Edge &other ) const
        {
            return _src == other._src && _trgt == other._trgt;
        }

        const Vertex &source() const
        {
            return _src;
        }

        const Vertex &target() const
        {
            return _trgt;
        }
    };

    struct GraphHashing
    {
        std::size_t operator()( const Vertex &v ) const
        {
            return std::hash< V >{}( v._data );
        }
        std::size_t operator()( const Edge &e ) const
        {
            auto h1 = this->operator ()( e._src );
            auto h2 = this->operator ()( e._trgt );
            return h1 ^ ( h2 << 1 );
        }
    };

    struct GraphComparators
    {
        constexpr bool operator()(const Edge &lhs, const Edge &rhs) const
        {
            return lhs._src._data + lhs._trgt._data <
                    rhs._src._data + rhs._trgt._data;
        }

        constexpr bool operator()(const Vertex &lhs, const Vertex &rhs) const
        {
            return lhs._data < rhs._data;
        }
    };

    using Edges = std::vector< Edge >;
    using Connections = std::pair< Edges , Edges >;
    using TraversableOnceGraph = std::map< const Vertex , std::list< Edge > , Graph< V >::GraphComparators>;
    using Path = std::list< Vertex >;

    void addEdges( const Edges &edges )
    {
        for( const auto &e : edges )
        {
            auto source = e._src;
            auto target = e._trgt;
            _graph[ source ].second.push_back( e );
            _graph[ target ].first .push_back( e );
        }
    }

    void addEdges( const Vertex &v , const Edges &edgesIn , const Edges &edgesOut )
    {
        auto &p = _graph[ v ] ;
        for( const auto &e : edgesIn )
            p.first.push_back( e );
        for( const auto &e : edgesOut )
            p.second.push_back( e );
    }

    const auto &vertex( const Vertex &v ) const
    {
        return _graph.at( v );
    }

    const auto &edgesIn( const Vertex &v ) const
    {
        return _graph.at( v ).first;
    }

    const auto &edgesOut( const Vertex &v ) const
    {
        return _graph.at( v ).second;
    }

    std::size_t edgesCount() const
    {
        return std::accumulate( _graph.begin() , _graph.end() , 0 ,
                                []( std::size_t count , const std::pair< Vertex , Connections > &p )
        {
            return count + p.second.first.size() ;
        });
    }

    TraversableOnceGraph getTraversableOnceGraph() const
    {
        TraversableOnceGraph graph;
        for( const auto &p : _graph )
            graph[ p.first ] = std::list< Edge >( p.second.second.begin() ,
                                                  p.second.second.end());
        return graph;
    }
    const auto &data() const
    {
        return _graph;
    }

    static std::list< Vertex > depthFirstTraverse( TraversableOnceGraph &graph , Vertex &seed )
    {
        Vertex currentVertex = seed;
        Path path;
        while(1)
        {
            std::list< Edge > &edgesOut = graph.at( currentVertex );
            if( edgesOut.empty())
                break;

            Edge nextEdge = edgesOut.front();
            edgesOut.pop_front();
            currentVertex = nextEdge.target();
            path.push_back( currentVertex );
        }
        return path;
    }

    static std::list< Vertex > extractEulerianCycle( TraversableOnceGraph &graph ,
                                                     const Vertex *seed = nullptr )
    {
        Path eCycle({ ((seed == nullptr) ? graph.begin()->first : *seed )});
        typename Path::iterator insertPosition = eCycle.begin();
        while( insertPosition != eCycle.end())
        {
            Path path = depthFirstTraverse( graph , *(insertPosition++) );
            eCycle.insert( insertPosition , path.begin() , path.end());
            insertPosition = std::find_if( eCycle.begin() , eCycle.end() ,
                                           [&graph]( const Vertex &v )
            {
                return !graph.at( v ).empty();
            });
        }
        return eCycle;
    }

    static bool emptyGraph( const TraversableOnceGraph &graph )
    {
        return std::all_of( graph.begin() , graph.end() ,
                            []( const std::pair< Vertex , Edges > &p )
        {
            return p.second.empty();
        });
    }

    static typename TraversableOnceGraph::iterator getVertexWithUnexpandedEdge( const TraversableOnceGraph &graph )
    {
        return std::find_if( graph.begin() , graph.end() ,
                             []( const std::pair< Vertex , Edges > &p )
        {
            return !p.second.empty();
        });
    }

    std::list< Vertex > extractEulerianPath() const
    {
        auto seedPoint = std::find_if( _graph.begin() , _graph.end() ,
                                       []( const std::pair< Vertex , Connections > &p )
        {
            // Out edges - In edges
            int diff = p.second.second.size() - p.second.first.size();
            return diff > 0 && diff % 2 == 1;
        });

        auto traversableOnceGraph = getTraversableOnceGraph();
        if( seedPoint == _graph.end())
            return extractEulerianCycle( traversableOnceGraph );
        else
            return extractEulerianCycle( traversableOnceGraph ,
                                         std::addressof( seedPoint->first ));
    }

    static std::string toString( const Path &path )
    {
        std::list< std::string > out;
        std::transform( std::begin( path ) , std::end( path ) ,
                        std::inserter( out , std::end( out )) ,
                        []( const rosalind::ba3::Graph< int >::Vertex &v )
        {
            return std::to_string( v.data());
        });
        return rosalind::io::join( out , "->");
    }


private:
    std::map< Vertex , Connections , GraphComparators > _graph;
};




/**
 * ba3b
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
 * ba3c
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
 * ba3d
 * @brief constructDeBruijnGraph
 * @param firstIt
 * @param lastIt
 * @return
 */
template< typename SeqIt >
std::map< const std::string , std::vector< std::string >>
constructDeBruijnGraph( SeqIt firstIt , SeqIt lastIt )
{
    // Construct the overlap map, O(n).
    auto prefixSuffixMap = constructOverlapMap( firstIt , lastIt );

    // Construct the graph, O(n).
    std::map< const std::string , std::vector< std::string >> deBruijnGraph;
    for( const auto &prefixSuffix : prefixSuffixMap )
        for( const auto prefix : prefixSuffix.second.first )
        {
            auto prefixStr = prefixSuffix.first;
            prefixStr.insert( 0 , 1 , prefix );
            prefixStr.pop_back();
            deBruijnGraph[ prefixStr ].emplace_back( prefixSuffix.first );
        }
    return deBruijnGraph;
}

/**
 * ba3e
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

/**
 * ba3f
 * @brief findEulerianPath
 * @param first
 * @param last
 */
template< typename SeqIt >
auto
findEulerianCycle( SeqIt first , SeqIt last )
{
    Graph< int > graph;
    using E = Graph< int >::Edge;
    using V = Graph< int >::Vertex;
    std::vector< E > edges;

    std::for_each( first , last ,
                   [&edges]( const std::string &line )
    {
        auto link = rosalind::io::split( line , " -> ");
        V src( atoi( link[0].c_str()));
        for( auto target : rosalind::io::split( link[1] , ','))
        {
            E e( src , V( atoi( target.c_str())));
            edges.push_back( e );
        }
    });
    graph.addEdges( edges );
    auto traversableOnceGraph = graph.getTraversableOnceGraph();
    return Graph< int >::extractEulerianCycle( traversableOnceGraph );
}

/**
 * ba3g
 * @brief findEulerianPath
 * @param first
 * @param last
 */
template< typename SeqIt >
auto
findEulerianPath( SeqIt first , SeqIt last )
{
    Graph< int > graph;
    using E = Graph< int >::Edge;
    using V = Graph< int >::Vertex;
    std::vector< E > edges;

    std::for_each( first , last ,
                   [&edges]( const std::string &line )
    {
        auto link = rosalind::io::split( line , " -> ");
        V src( atoi( link[0].c_str()));
        for( auto target : rosalind::io::split( link[1] , ','))
        {
            E e( src , V( atoi( target.c_str())));
            edges.push_back( e );
        }
    });
    graph.addEdges( edges );
    return graph.extractEulerianPath();
}


}
}
#endif // BA2_HPP
