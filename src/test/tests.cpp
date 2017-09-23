#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()

#include "tests.h"
#include "catch/catch.hpp"


TEST_CASE("Basic Utilities")
{
    SECTION("JOIN: joining string items with separator.")
    {
        REQUIRE( rosalind::io::join({"item1","item2","item3"},",") ==
                 "item1,item2,item3");
        REQUIRE( rosalind::io::join({"apple","banana"},"-") ==
                 "apple-banana");
        REQUIRE( rosalind::io::join({"hello"} , "&") == "hello");

        REQUIRE( rosalind::io::join({},",") == "");
    }

    SECTION("Set Based Equality: check the contents "
            "equality of two containers regardless the order of elements.")
    {
        using namespace test_utils;
        using V = std::vector< int >;
        REQUIRE( setBasedEquality( V{1,3,2} , V{3,1,2} ));
        REQUIRE( setBasedEquality( V{3,3,3,1} , V{1,3,1} , false ));
        REQUIRE(!setBasedEquality( V{3,2} , V{3,4} ));
    }

    SECTION("Split String with string delimeter")
    {
        using S = std::string;
        using V = std::vector< S >;
        REQUIRE( rosalind::io::split( "AACCCA#-#GTTTGA" , "#-#") ==
                 V({"AACCCA","GTTTGA"}) );

        REQUIRE( rosalind::io::split( "AACCCA -> GTTTGA" , " -> ") ==
                 V({"AACCCA","GTTTGA"}) );
    }

    SECTION("Find with mismatches")
    {
        //TODO
    }
}


TEST_CASE("Finding Hidden Messages Algorithms","[BA1]")
{
    using namespace test_utils;
    SECTION("BA1A: Pattern Count")
    {
        auto inputLines = getFileLines( dataFilePath("ba1a"));
        REQUIRE( rosalind::ba1::patternCount( inputLines[ 0 ] , inputLines[ 1 ]) ==
                atoi( getFileLines( outputFilePath( "ba1a" ))[0].c_str( )));
    }

    SECTION("BA1B: Most Frequent Words")
    {
        auto inputLines = getFileLines( dataFilePath("ba1b"));

        auto actual =
                rosalind::ba1::frequentWordsBruteForce( inputLines );
        auto expected =
                rosalind::io::split( getFileLines(
                                         outputFilePath( "ba1b" ))[0] , ' ');

        REQUIRE( expected.size() == actual.size());
        REQUIRE( actual.size() != 0 );
        REQUIRE( setBasedEquality( actual , expected ));
    }

    SECTION("BA1C: Complement Sequence")
    {
        auto inputLines = getFileLines( dataFilePath("ba1c"));
        auto expected = getFileLines(outputFilePath( "ba1c" ))[0];
        auto actual = rosalind::ba1::complementSequence( inputLines[0] );
        auto actual2 = rosalind::ba1::complementSequence2( inputLines[0] );
        REQUIRE( expected == actual );
        REQUIRE( expected == actual2 );
    }

    SECTION("BA1D: Pattern Matching")
    {
        auto inputLines = getFileLines( dataFilePath("ba1d"));
        auto _expected =
                rosalind::io::split( getFileLines(
                                         outputFilePath( "ba1d" ))[0] , ' ');

        auto actual = rosalind::ba1::patternMatching( inputLines );
        decltype( actual ) expected;
        std::transform( _expected.begin() , _expected.end() ,
                        std::inserter( expected , expected.begin()) ,
                        []( const std::string &s ) { return std::atoi( s.c_str()); } );
        REQUIRE( expected == actual );
    }

    SECTION( "BA1E: Finding Clumps Problem" )
    {
        auto inputLines = getFileLines( dataFilePath( "ba1e"));
        REQUIRE( inputLines.size() == 2 );

        auto inputSequence = inputLines[0];
        auto parametersStr = rosalind::io::split( inputLines[1] , ' ');
        REQUIRE( parametersStr.size() == 3 );

        auto expectedOutput =
                rosalind::io::split( getFileLines(
                                         outputFilePath( "ba1e" ))[0] , ' ');
        auto actualOutput   =
                rosalind::ba1::findClumps( inputLines );

        REQUIRE( setBasedEquality( expectedOutput , actualOutput ));
    }

    SECTION( "BA1F: Find a Position in a Genome Minimizing the Skew")
    {
        auto inputLines =
                getFileLines( test_utils::dataFilePath("ba1f"));
        auto expectedOutput =
                getFileLines( test_utils::outputFilePath("ba1f"));
        auto actualOutput =
                rosalind::ba1::skewDiagram( inputLines[0] );

        auto _expected = rosalind::io::split( expectedOutput[0] , ' ');
        std::set< int > expected;
        std::transform( _expected.begin() , _expected.end() ,
                        std::inserter( expected , expected.begin()) ,
                        []( const std::string &s ) { return std::atoi( s.c_str()); } );

        auto _ = actualOutput.first.first;
        auto actual = std::set< int >( _.begin() , _.end());

        REQUIRE( actual == expected );
    }

    SECTION( "BA1G: Hamming Distance" )
    {
        auto input = getFileLines( test_utils::dataFilePath("ba1g"));
        auto expected = atoi( getFileLines( test_utils::outputFilePath("ba1g"))[0].c_str());
        using namespace rosalind::ba1;
        REQUIRE( hammingDistance( input[0] , input[1] ) == expected );
        REQUIRE( hammingDistance( input[0].c_str() , input[1].c_str() , input[0].size() ) ==
                expected );
    }

    SECTION( "BA1H: Approximate Pattern Matching" )
    {
        auto input = getFileLines( dataFilePath("ba1h"));
        auto actual = rosalind::ba1::approximatePatternMatching( input );

        auto _expected =
                rosalind::io::split( getFileLines(
                                         outputFilePath( "ba1h" ))[0] , ' ');

        decltype( actual ) expected;
        std::transform( _expected.begin() , _expected.end() ,
                        std::inserter( expected , expected.begin()) ,
                        []( const std::string &s ) { return std::atoi( s.c_str()); } );
        REQUIRE( expected == actual );
    }

    SECTION( "BA1I: Most Frequent Words With Mismatches")
    {
        auto input = getFileLines( dataFilePath("ba1i"));
        auto actual = rosalind::ba1::frequentWordsWithMismatches( input );
        auto expected =
                rosalind::io::split( getFileLines(
                                         outputFilePath( "ba1i" ))[0] , ' ');
        REQUIRE( setBasedEquality( actual , expected ));
    }

    SECTION( "BA1J: Most Frequent Words With Mismatches And Reverse Complement" )
    {
        auto input = getFileLines( dataFilePath("ba1j"));
        auto actual = rosalind::ba1::frequentWordsWithMismatchesAndReverseComplement( input );
        auto expected =
                rosalind::io::split( getFileLines(
                                         outputFilePath( "ba1j" ))[0] , ' ');
        REQUIRE( setBasedEquality( actual , expected ));
    }

    SECTION( "BA1K: String Frequency Array")
    {
        auto input = getFileLines( dataFilePath("ba1k"));
        auto actual = rosalind::ba1::stringFrequencyArray( input );
        auto _expected =
                rosalind::io::split( getFileLines(
                                         outputFilePath( "ba1k" ))[0] , ' ');

        decltype( actual ) expected;
        std::transform( _expected.begin() , _expected.end() ,
                        std::inserter( expected , expected.begin()) ,
                        []( const std::string &s ) { return std::atoi( s.c_str()); } );
        REQUIRE( expected == actual );
    }

    SECTION( "BA1L")
    {
        auto input = getFileLines( dataFilePath("ba1l"));
        auto actual = rosalind::ba1::encode( input );
        auto expected = std::atoll( getFileLines( outputFilePath( "ba1l" ))[0].c_str());
        REQUIRE( actual == expected );
    }

    SECTION( "BA1M")
    {
        auto input = getFileLines( dataFilePath("ba1m"));
        auto actual = rosalind::ba1::numberToPattern( input );
        auto expected = getFileLines( outputFilePath( "ba1m" ))[0];
        REQUIRE( actual == expected );
    }

    SECTION( "BA1N")
    {
        auto input = getFileLines( dataFilePath("ba1n"));
        auto actual = rosalind::ba1::dNeighborhood( input );
        auto expected = getFileLines( outputFilePath("ba1n"));
        REQUIRE( setBasedEquality( actual , expected ));
    }


    SECTION( "BA5A: Find the minimum number of coins needed to make change.")
    {
        auto input =
                getFileLines( dataFilePath("ba5a"));
        auto _expectedOutput =
                getFileLines( outputFilePath("ba5a"));
        auto expectedOutput =
                std::atoi( _expectedOutput[ 0 ].c_str());

        REQUIRE( rosalind::ba1::minimumCoinsChange( input ) ==
                 expectedOutput );
    }

}


TEST_CASE("Finding Motifs Algorithms","[BA2]")
{
    using namespace test_utils;
    SECTION("BA2A: Motif Enumeration")
    {
        auto input = getFileLines( dataFilePath("ba2a"));
        auto actual = rosalind::ba2::motifEnumeration( input );
        auto expected = getFileLines( outputFilePath("ba2a"));
        REQUIRE( setBasedEquality( expected , actual ));
    }

    SECTION("BA2B: Median Kmer")
    {
        auto input = getFileLines( dataFilePath("ba2b"));
        auto actual = rosalind::ba2::medianKmer( input );
        auto expected = getFileLines( outputFilePath("ba2b"));
        REQUIRE( expected[0] == actual );
    }

    SECTION("BA2C: Profile-most probable kmer")
    {
        auto input = getFileLines( dataFilePath("ba2c"));
        auto actual = rosalind::ba2::profileMostProbableKmer( input );
        auto expected = getFileLines( outputFilePath("ba2c"));
        REQUIRE( expected[0] == actual );
    }

    SECTION("BA2D: Greedy Motif Search")
    {
        auto input = getFileLines( dataFilePath("ba2d"));
        auto actual = rosalind::ba2::greedyMotifSearch( input );
        auto expected = getFileLines( outputFilePath("ba2d"));
        REQUIRE( expected == actual );
    }

    SECTION("BA2E: Greedy Motif Search With Pseudocounts")
    {
        auto input = getFileLines( dataFilePath("ba2e"));
        auto actual = rosalind::ba2::greedyMotifSearchWithPseudoCount( input );
        auto expected = getFileLines( outputFilePath("ba2e"));
        REQUIRE( expected == actual );
    }

    SECTION("BA2F: Randomized Motif Search")
    {
        auto input = getFileLines( dataFilePath("ba2f"));
        auto actual = rosalind::ba2::randomizedMotifSearch( input );
        auto expected = getFileLines( outputFilePath("ba2f"));
        REQUIRE( expected == actual );
    }

    SECTION("BA2G: Gibbs Sampler Motif Finding")
    {
        auto input = getFileLines( dataFilePath("ba2g"));
        auto actual = rosalind::ba2::gibbsSampler( input );
        auto expected = getFileLines( outputFilePath("ba2g"));

        CAPTURE( rosalind::ba2::hammingDistanceScore( actual ));
        CAPTURE( rosalind::ba2::hammingDistanceScore( expected ));

        REQUIRE( expected == actual );
    }

    SECTION("BA2H: Distance Between Pattern And Strings")
    {
        auto input = getFileLines( dataFilePath("ba2h"));
        auto actual = rosalind::ba2::minimalHammingDistance( input );
        auto expected = atoi( getFileLines( outputFilePath("ba2h"))[0].c_str());
        REQUIRE( expected == actual );
    }
}

TEST_CASE("Graph Algorithms","[BA3]")
{
    using namespace test_utils;
    SECTION("BA3a: Find Kmer Composition of String")
    {
        auto input = getFileLines( dataFilePath("ba3a"));
        auto actual = rosalind::ba1::extractKmers( input[1] , atoi( input[0].c_str()));
        auto expected = getFileLines( outputFilePath("ba3a"));
        REQUIRE( setBasedEquality( actual , expected ));
    }

    SECTION("BA3b: String Reconstruction From Kmers")
    {
        auto input = getFileLines( dataFilePath("ba3b"));
        auto actual = rosalind::ba3::reconstructStringFromKmers(
                    input.cbegin() , input.cend());
        auto expected = getFileLines( outputFilePath("ba3b"))[0];
        REQUIRE( actual == expected );
    }

    SECTION("BA3c: Overlap Strings Graph")
    {
        auto input = getFileLines( dataFilePath("ba3c"));
        auto actual = rosalind::ba3::overlapStringsGraph( input.cbegin() , input.cend());
        auto expected = getFileLines( outputFilePath("ba3c"));
        decltype( expected ) _actual;
        std::transform( std::begin( actual ) , std::end( actual ) ,
                        std::inserter( _actual , _actual.begin()) ,
                        []( const std::pair< std::string , std::string > &p )
        {
            return p.first + " -> " + p.second;
        });
        REQUIRE( setBasedEquality( _actual , expected ));
    }

    SECTION("BA3d: De Bruijn Overlap Strings Graph")
    {
        using P = std::pair< std::string , std::vector< std::string >>;
        auto input = getFileLines( dataFilePath("ba3d"));
        auto actual = rosalind::ba3::constructDeBruijnGraph( input );

        auto expected = getFileLines( outputFilePath("ba3d"));
        decltype( expected ) _actual;
        std::transform( std::begin( actual ) , std::end( actual ) ,
                        std::inserter( _actual , _actual.end()) ,
                        []( const std::pair< std::string , std::vector< std::string >> &p )
        {
            auto &targets = p.second;
            return p.first + " -> " + rosalind::io::join( targets , ",");
        });

        CAPTURE( expected.size());
        CAPTURE( _actual.size());

        std::set< std::string > expectedSet( std::begin( expected ) , std::end( expected ));
        std::set< std::string > actualSet( std::begin( _actual ) , std::end( _actual ));
        auto firstMismatch = [expected,_actual](){
            auto p = setBasedFirstMismatch( expected , _actual );
            return std::string("\nExpected:") + p.first + "\nActual:" + p.second; };

        CAPTURE( firstMismatch());
//        REQUIRE( setBasedEquality( expected , _actual ));
    }

    SECTION("BA3e: De Bruijn Overlap Strings Graph")
    {
        using S = std::string;
        using V = std::vector< S >;
        using P = std::pair< const S , std::vector< S >>;
        using M = std::map< const S , std::vector< S >>;

        auto input = getFileLines( dataFilePath("ba3e"));
        auto actual = rosalind::ba3::constructDeBruijnGraph( input.cbegin() , input.cend());
        V _actual;
        std::transform( std::begin( actual ) , std::end( actual ) ,
                        std::inserter( _actual , std::end( _actual )) ,
                        []( const P &p )
        {
            auto _p = p;
            std::sort( _p.second.begin() , _p.second.end());
            return _p.first + " -> " + rosalind::io::join( _p.second , ",");
        });

        auto expected = getFileLines( outputFilePath("ba3e"));
        decltype( expected ) _expected;
        std::transform( std::begin( expected ) , std::end( expected ) ,
                        std::inserter( _expected , _expected.end()) ,
                        []( const S &s )
        {
            auto l = rosalind::io::split( s , " -> ");
            auto k = l[0];
            auto v = rosalind::io::split(l[1] , ',');
            std::sort( v.begin() , v.end());
            return k + " -> " + rosalind::io::join( v , ",");
        });

        CAPTURE( _expected.size());
        CAPTURE( _actual.size());
        REQUIRE( setBasedEquality( _expected , _actual )) ;
    }

    SECTION("BA3f: Find Eulerian Cycle")
    {
        /**
          * @todo correctness check method.
          **/

    }

    SECTION("BA3g: Find Eulerian Path")
    {
        /**
          * @todo correctness check method.
          **/

    }

    SECTION("BA3h: String Reconstruction From Kmers")
    {
        /**
          * @todo correctness check method.
          **/

    }
    SECTION("BA3i: ")
    {
        /**
          * @todo correctness check method.
          **/

    }
    SECTION("BA3j: ")
    {
        auto input = getFileLines( dataFilePath("ba3j"));
        unsigned int d = atoi( rosalind::io::split( input.front() , " ")[1].c_str() );
        auto actual = rosalind::ba3::reconstructStringFromSparsePairedKmers(
                    input.cbegin() + 1 , input.cend(), d );
        auto expected = getFileLines( outputFilePath("ba3j")).front();
        REQUIRE( actual == expected );
    }
    SECTION("BA3jk: Contigs generation")
    {
        auto input = getFileLines( dataFilePath("ba3k"));
        auto actual = rosalind::ba3::generateContigs( input.cbegin() , input.cend());
        auto expected = rosalind::io::split(
                    getFileLines( outputFilePath("ba3k")).front() , " ");
        REQUIRE( setBasedEquality( actual , expected ));
    }
    SECTION("BA3l: Construct a String Spelled by a Gapped Genome Path")
    {
        auto input = getFileLines( dataFilePath("ba3l"));
        auto k = std::stoi( rosalind::io::split( input.front() , " ").front());
        auto d = std::stoi( rosalind::io::split( input.front() , " ").at( 1 ));
        auto actual = rosalind::ba3::constructPairedStringsFromGappedPath(
                    input.cbegin() + 1 , input.cend() , k , d );
        auto expected = getFileLines( outputFilePath("ba3l")).front();
        REQUIRE( actual == expected );
    }
    /**
    SECTION("BA3m: Generate All Maximal Non-Branching Paths in a Graph")
    {
        auto input = getFileLines( dataFilePath("ba3m"));
        auto actual = rosalind::ba3::getMaximalNonBranchingPaths(
                    input.cbegin(), input.cend());
        auto expected = getFileLines( outputFilePath("ba3m"));
        REQUIRE( setBasedEquality( actual , expected ));
    }**/
}

TEST_CASE("Protein Processing & Analysis","[BA3]")
{
    using namespace test_utils;
    SECTION("BA4a: Translate an RNA String into an Amino Acid String")
    {
        auto input = getFileLines( dataFilePath("ba4a"));
        auto actual = rosalind::ba4::translateRNA2AA(
                    input.front().cbegin() , input.front().cend() );
        auto expected = getFileLines( outputFilePath("ba4a")).front();
        REQUIRE( actual == expected );
    }
    SECTION("BA4b: Find Substrings of a Genome Encoding a Given Amino Acid String")
    {
        auto input = getFileLines( dataFilePath("ba4b"));
        auto actual = rosalind::ba4::substringEncodingAA(
                    input.front().cbegin() , input.front().cend() ,
                    input.at( 1 ));
        auto expected = getFileLines( outputFilePath("ba4b"));
        REQUIRE( setBasedEquality( actual , expected ));
    }
}
