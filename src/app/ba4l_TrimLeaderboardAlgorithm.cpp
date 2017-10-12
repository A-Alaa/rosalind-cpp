#include "ba4.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    std::vector< std::string > input;
    if( argc > 1 )
        input = rosalind::io::getFileLines( argv[1] );
    else input = rosalind::io::readInputStream();

    int n = std::stoi( input.at( 2 ));
    std::vector< int > spectrum;
    for( const std::string &strMass : rosalind::io::split( input.at( 1 ) , " "))
        spectrum.push_back( std::stoi( strMass ));

    std::multimap< int , std::pair< std::string , int > , std::greater< int >> leaderboard;

    auto score = [spectrum]( const std::vector< uint8_t > &p )->int{
        return rosalind::ba4::scoreLinearPeptideSpectrum( p , spectrum.cbegin() , spectrum.cend() );
    };
    for( std::string &peptide : rosalind::io::split( input.front() , " "))
    {
        std::vector< uint8_t > pept;
        std::transform( peptide.cbegin() , peptide.cend() ,
                        std::back_inserter( pept ) ,
                        []( char aa ){
            return rosalind::ba4::massTable.at( aa );
        });
        int scorePept = score( pept );
        leaderboard.emplace( scorePept , std::make_pair( std::move( peptide ) , 0 ));
    }

    rosalind::ba4::trimLeaderboard( leaderboard , n , []( int m ){ return true;} );
    std::vector< std::string > leadingPeptides;
    std::transform( leaderboard.cbegin() , leaderboard.cend() ,
                    std::back_inserter( leadingPeptides ) ,
                    []( const std::pair< int , std::pair< std::string , int > > &p ){
        std::cout << p.first << std::endl;
        return p.second.first;
    });

    std::cout << rosalind::io::join( leadingPeptides.cbegin() , leadingPeptides.cend() , " " )
              << std::endl;
    return 0;
}
