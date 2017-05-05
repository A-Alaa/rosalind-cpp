#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()

#include "tests.h"
#include "catch/catch.hpp"

TEST_CASE( "BA1E: Finding Clumps Problem", "[BA1E]" ) {
    auto inputLines = test_utils::getFileLines( test_utils::dataFilePath( "ba1e"));
    REQUIRE( inputLines.size() == 2 );

    auto inputSequence = inputLines[0];
    auto parametersStr = test_utils::split( inputLines[1] , ' ');
    REQUIRE( parametersStr.size() == 3 );

    auto expectedOutput =
            test_utils::split( test_utils::getFileLines(
                                   test_utils::outputFilePath( "ba1e" ))[0] , ' ');
    auto actualOutput   =
            rosalind::basic::findClumps(    inputSequence ,
                                            std::atoi( parametersStr[0].c_str( )) ,
                                            std::atoi( parametersStr[1].c_str( )) ,
                                            std::atoi( parametersStr[2].c_str( )));

    const auto expected = std::set< std::string >( expectedOutput.begin() , expectedOutput.end());
    const auto actual    = std::set< std::string >( actualOutput.begin()   , actualOutput.end());
    REQUIRE( expected == actual );
}


TEST_CASE( "BA1F: Find a Position in a Genome Minimizing the Skew" , "[BA1F]")
{
    auto inputLines = test_utils::getFileLines( test_utils::dataFilePath("ba1f"));
    auto expectedOutput = test_utils::getFileLines( test_utils::outputFilePath("ba1f"));
    auto actualOutput =
            rosalind::basic::skewDiagram( inputLines[0] );

    auto _expected = test_utils::split( expectedOutput[0] , ' ');
    std::set< int > expected;
    std::transform( _expected.begin() , _expected.end() ,
                    std::inserter( expected , expected.begin()) ,
                    []( const std::string &s ) { return std::atoi( s.c_str()); } );

    auto _ = actualOutput.first.first;
    auto actual = std::set< int >( _.begin() , _.end());

    REQUIRE( actual == expected );
}

