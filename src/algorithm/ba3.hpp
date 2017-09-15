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

        operator V() const
        {
            return _data;
        }

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
        Edge( const V &srcData , const V &trgtData )
            : _src( srcData ) , _trgt( trgtData ) {}

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
        bool operator()(const Edge &lhs, const Edge &rhs) const
        {
            return lhs._src._data + lhs._trgt._data <
                    rhs._src._data + rhs._trgt._data;
        }

        bool operator()(const Vertex &lhs, const Vertex &rhs) const
        {
            return lhs._data < rhs._data;
        }
    };

    using Edges = std::vector< Edge >;
    using Connections = std::pair< Edges , Edges >;
    using TraversableOnceGraph = std::map< const Vertex , std::list< Edge > , GraphComparators>;
    using Path = std::list< Vertex >;

    Graph() = default;

    Graph( const Edges &edges )
    {
        addEdges( edges );
    }

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

    std::string toString() const
    {
        std::string str;
        for( const auto &p : _graph )
        {
            str += p.first;
            str += "->";
            std::vector< std::string > targets;
            std::transform( p.second.second.begin() ,
                            p.second.second.end() ,
                            std::inserter( targets , targets.end()) ,
                            []( const Edge &e ){
                return  e.target().data();
            });
            str += io::join( targets , "," );
            str += '\n';
        }
        return str;
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
    str.reserve( firstIt->length() + std::distance( std::next( firstIt , 1 ) , lastIt ));
    str.insert( str.begin() , firstIt->begin() , firstIt->end());
    std::transform( std::next( firstIt , 1 ) , lastIt ,
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

/**
 * ba3h
 * @brief reconstructStringFromSparseKmers
 * @param first
 * @param last
 * @return
 */
template< typename SeqIt >
std::string
reconstructStringFromSparseKmers( SeqIt first , SeqIt last )
{
    using G = Graph< std::string >;
    using E = G::Edge;
    using V = G::Vertex ;

    auto deBruijnGraph = constructDeBruijnGraph( first , last );
    G::Edges edges;
    std::for_each( deBruijnGraph.begin() , deBruijnGraph.end() ,
                   [&edges]( const std::pair< std::string , std::vector< std::string >> &p )
    {
        V src( p.first );
        for( const auto &target : p.second )
        {
            E e( src , V( target ));
            edges.push_back( e );
        }
    });
    deBruijnGraph.clear();

    auto graph = G( edges );
    edges.clear();

    auto eulerianPath = graph.extractEulerianPath();
    std::vector< std::string > kmers;
    std::transform( eulerianPath.begin() , eulerianPath.end() ,
                    std::inserter( kmers , kmers.end()) ,
                    []( const V &v )
    {
        return v.data();
    });
    eulerianPath.clear();

    return reconstructStringFromKmers( kmers.begin() , kmers.end());

}

template< typename SeqIt >
auto
constructPairedStringsOverlapMap( SeqIt firstIt , SeqIt lastIt )
{
    using S = std::string;
    using V = std::vector< std::array< char , 2 >>;
    std::unordered_map< S , std::pair< V , V >> prefixSuffixMap;
    const unsigned int k = firstIt->size() / 2;
    // Populate suffixes, O(n).
    for( auto it = firstIt ; it != lastIt ; ++it )
    {
        // Remove chars at x in string: xCCC|xCCC
        S suffix = S( it->begin() + 1 , it->begin() + k + 1 ) +
                S( it->begin() + k + 2 , it->end());
        assert( suffix.size() == 2 * k - 1 );
        prefixSuffixMap[ suffix ].first.push_back(
        {it->front(), it->at( k + 1 )});
    }

    // Populate prefixes, O(n).
    for( auto it = firstIt ; it != lastIt ; ++it )
    {
        // Remove chars at x in string: CCCx|CCCx
        S prefix = S( it->begin() , it->begin() + k - 1 ) +
                S( it->begin() + k , it->end() - 1 );
        assert( prefix.size() == 2 * k - 1 );
        prefixSuffixMap[ prefix ].second.push_back(
        { it->at( k - 1 ) , it->back() });
    }

    return prefixSuffixMap;
}

template< typename SeqIt >
std::map< const std::string , std::vector< std::string >>
constructPairedDeBruijnGraph( SeqIt firstIt , SeqIt lastIt , char sep = '|' )
{
    auto checkInput = [firstIt,lastIt,sep](){
        unsigned int k = firstIt->size() / 2;
        return std::all_of( firstIt , lastIt ,
                            [k,sep]( const std::string &p ){
            return p.size() == 2 * k + 1 &&
                    p[k] == sep;
        });
    };
    assert( checkInput());
    const unsigned int k = firstIt->size() / 2;

    // Construct the overlap map, O(n).
    auto prefixSuffixMap = constructPairedStringsOverlapMap( firstIt , lastIt );

    // Construct the graph, O(n).
    std::map< const std::string , std::vector< std::string >> deBruijnGraph;
    for( const auto &prefixSuffix : prefixSuffixMap )
        for( const std::array< char , 2 > &prefix : prefixSuffix.second.first )
        {
            // prefixStr = CCx|CCx
            auto prefixStr = prefixSuffix.first;
            // add prefixes so prefixStr = pCCx|pCCx
            prefixStr.insert( 0     , 1 , prefix[0] );
            prefixStr.insert( k + 1 , 1 , prefix[1] );
            // now remove x's
            prefixStr.pop_back();
            prefixStr.erase( k - 1 , 1 );
            deBruijnGraph[ prefixStr ].emplace_back( prefixSuffix.first );
        }
    return deBruijnGraph;
}


template< typename SeqIt >
std::string
reconstructStringFromSparsePairedKmers( SeqIt first , SeqIt last ,
                                        unsigned int spaced )
{
    using G = Graph< std::string >;
    using E = G::Edge;
    using V = G::Vertex ;
    const int k = io::split( *first , '|').front().size();
    auto deBruijnGraph = constructPairedDeBruijnGraph( first , last );
    G::Edges edges;
    std::for_each( deBruijnGraph.begin() , deBruijnGraph.end() ,
                   [&edges]( const std::pair< std::string , std::vector< std::string >> &p )
    {
        V src( p.first );
        for( const auto &target : p.second )
        {
            E e( src , V( target ));
            edges.push_back( e );
        }
    });
    deBruijnGraph.clear();

    auto graph = G( edges );
    edges.clear();

    auto eulerianPath = graph.extractEulerianPath();
    std::vector< std::string > kmers;
    std::transform( eulerianPath.cbegin() ,
                    eulerianPath.cend() ,
                    std::inserter( kmers , kmers.end()) ,
                    []( const V &v )
    {
        return io::split( v.data() , "|" ).front();
    });
    // Append the residual part.
    std::transform( std::prev( eulerianPath.cend() , spaced + k ) ,
                    eulerianPath.cend() ,
                    std::inserter( kmers , kmers.end()) ,
                    []( const V &v )
    {
        return io::split( v.data() , "|" ).at( 1 );
    });
    return reconstructStringFromKmers( kmers.cbegin() , kmers.cend());

}
}
}
#endif // BA2_HPP
