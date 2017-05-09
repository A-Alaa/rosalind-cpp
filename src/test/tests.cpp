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

