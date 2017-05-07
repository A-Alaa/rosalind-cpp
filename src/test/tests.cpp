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
}

TEST_CASE("Finding Hidden Messages Algorithms","[BA1]")
{
    using namespace test_utils;
    SECTION("BA1A: Pattern Count")
    {
        auto inputLines = getFileLines( dataFilePath("ba1a"));
        REQUIRE( rosalind::basic::patternCount( inputLines[ 0 ] , inputLines[ 1 ]) ==
                atoi( getFileLines( outputFilePath( "ba1a" ))[0].c_str( )));
    }

    SECTION("BA1B: Most Frequent Words")
    {
        auto inputLines = getFileLines( dataFilePath("ba1b"));

        auto actualOutput =
                rosalind::basic::frequentWordsBruteForce( inputLines );
        auto expectedOutput =
                rosalind::io::split( getFileLines(
                                         outputFilePath( "ba1b" ))[0] , ' ');
        const auto expected = std::set< std::string >( expectedOutput.begin() ,
                                                       expectedOutput.end());
        const auto actual   = std::set< std::string >( actualOutput.begin() ,
                                                       actualOutput.end());

        REQUIRE( expected.size() == actual.size());
        REQUIRE( actual.size() != 0 );
        REQUIRE( expected == actual );
    }

    SECTION("BA1C: Complement Sequence")
    {
        auto inputLines = getFileLines( dataFilePath("ba1c"));
    }

    SECTION("BA1D: Pattern Matching")
    {
        auto inputLines = getFileLines( dataFilePath("ba1d"));
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
                rosalind::basic::findClumps( inputSequence ,
                                             std::atoi( parametersStr[0].c_str( )) ,
                std::atoi( parametersStr[1].c_str( )) ,
                std::atoi( parametersStr[2].c_str( )));

        const auto expected = std::set< std::string >( expectedOutput.begin() ,
                                                       expectedOutput.end());
        const auto actual   = std::set< std::string >( actualOutput.begin() ,
                                                       actualOutput.end());
        REQUIRE( expected == actual );
    }

    SECTION( "BA1F: Find a Position in a Genome Minimizing the Skew")
    {
        auto inputLines =
                getFileLines( test_utils::dataFilePath("ba1f"));
        auto expectedOutput =
                getFileLines( test_utils::outputFilePath("ba1f"));
        auto actualOutput =
                rosalind::basic::skewDiagram( inputLines[0] );

        auto _expected = rosalind::io::split( expectedOutput[0] , ' ');
        std::set< int > expected;
        std::transform( _expected.begin() , _expected.end() ,
                        std::inserter( expected , expected.begin()) ,
                        []( const std::string &s ) { return std::atoi( s.c_str()); } );

        auto _ = actualOutput.first.first;
        auto actual = std::set< int >( _.begin() , _.end());

        REQUIRE( actual == expected );
    }

    SECTION( "BA5A: Find the minimum number of coins needed to make change.")
    {
        auto input =
                getFileLines( dataFilePath("ba5a"));
        auto _expectedOutput =
                getFileLines( outputFilePath("ba5a"));
        auto expectedOutput =
                std::atoi( _expectedOutput[ 0 ].c_str());

        REQUIRE( rosalind::basic::minimumCoinsChange( input ) ==
                 expectedOutput );
    }

}

