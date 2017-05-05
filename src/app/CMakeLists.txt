## Dependencies
include( Qt )
include_directories( ../common )

## ba1e: FindingClumps
set( FINDING_CLUMPS_EXE "ba1e_FindingClumps" )
set( FINDING_CLUMPS_SOURCE_FILES
    ba1e_FindingClumps.cpp)

##  ba1f: MinimumSkew
set( MINIMUM_SKEW_EXE "ba1f_MinimumSkew" )
set( MINIMUM_SKEW_SOURCE_FILES
    ba1f_MinimumSkew.cpp)

## ba5a: MinCoinsChange
set( MINIMUM_COINS_CHANGE_EXE "ba5a_MinCoinsChange")
set( MINIMUM_COINS_CHANGE_SOURCE_FILES
    ba5a_MinCoinsChange.cpp)


### Generate the executable
add_executable( ${FINDING_CLUMPS_EXE} ${FINDING_CLUMPS_SOURCE_FILES})
add_executable( ${MINIMUM_SKEW_EXE} ${MINIMUM_SKEW_SOURCE_FILES})
add_executable( ${MINIMUM_COINS_CHANGE_EXE} ${MINIMUM_COINS_CHANGE_SOURCE_FILES})


# Install the executable (generate with 'make install' in the bin directory)
install( TARGETS ${FINDING_CLUMPS_EXE} DESTINATION bin )
