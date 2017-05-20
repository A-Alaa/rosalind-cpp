#ifndef TESTS_H
#define TESTS_H

#include "ba1.hpp"
#include "ba2.hpp"
#include "ba3.hpp"

#define DATA_DIRECTORY "/home/asem/GP/rosalind-cpp/data"
#define OUTPUT_DIRECTORY "/home/asem/GP/rosalind-cpp/output"


namespace test_utils
{
    auto getFileLines( const std::string &filePath )
    {
        std::ifstream f( filePath );
        std::vector< std::string > lines;
        std::string line;
        if( f )
            while( std::getline( f , line ))
                lines.push_back( line );
        else std::cout << "Failed to open file:" << filePath ;
        return lines;
    }

    auto constructFilePath( const std::string &directoryPath )
    {
        return [&]( const std::string &filename )
        {
            return directoryPath + "/rosalind_" + filename + ".txt";
        };
    }

    auto dataFilePath( const std::string &problemTag )
    {
        return constructFilePath( DATA_DIRECTORY )( problemTag );
    }

    auto outputFilePath( const std::string &problemTag )
    {
        return constructFilePath( OUTPUT_DIRECTORY )( problemTag );
    }

    template< typename Container1 , typename Container2 ,
              typename Elem = typename Container1::value_type >
    bool setBasedEquality( const Container1 &a , const Container2 &b ,
                           bool checkSize = true )
    {
        if( checkSize && a.size() != b.size())
            return false;

        std::set< Elem > setA( std::begin( a ) , std::end( a ));
        std::set< Elem > setB( std::begin( b ) , std::end( b ));
        return setA == setB;
    }

}
#endif // TESTS_H
