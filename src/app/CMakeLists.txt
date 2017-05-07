## Dependencies
#include( Qt )
include_directories( ../common )

### Generate the executable
add_executable( ba1a_PatternCount                       ba1a_PatternCount.cpp                       )
add_executable( ba1b_MostFrequentWords                  ba1b_MostFrequentWords.cpp                  )
add_executable( ba1c_ComplementSequence                 ba1c_ComplementSequence.cpp                 )
add_executable( ba1d_PatternMatching                    ba1d_PatternMatching.cpp                    )
add_executable( ba1e_FindingClumps                      ba1e_FindingClumps.cpp                      )
add_executable( ba1f_MinimumSkew                        ba1f_MinimumSkew.cpp                        )
add_executable( ba1g_HammingDistance                    ba1g_HammingDistance.cpp                    )
add_executable( ba1h_ApproximatePatternMatching         ba1h_ApproximatePatternMatching.cpp         )
add_executable( ba1i_MostFrequentWordsWithMismatches    ba1i_MostFrequentWordsWithMismatches.cpp    )

add_executable( ba5a_MinCoinsChange                     ba5a_MinCoinsChange.cpp                     )





# Install the executable (generate with 'make install' in the bin directory)
#install( TARGETS ${FINDING_CLUMPS_EXE} DESTINATION bin )
