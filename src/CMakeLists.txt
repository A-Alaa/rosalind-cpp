### Library name: lasso
set( ROSALIND_LIB "rosalind" )


### Library sources/headers
set( ROSALIND_LIB_SOURCE_FILES  )

set( ROSALIND_LIB_HEADER_FILES
    common/Common.hh
    common/catch.hpp )


### Include directories.
include_directories( . )
include_directories( common )
include_directories( catch )


#### Define the library
#add_library( ${ROSALIND_LIB}
#    ${ROSALIND_LIB_SOURCE_FILES}
#    ${ROSALIND_LIB_HEADER_FILES})

#set_target_properties( ${ROSALIND_LIB} PROPERTIES LINKER_LANGUAGE CXX COMPILE_DEFINITIONS "LIBRARY_API=1")

#install( TARGETS ${ROSALIND_LIB} DESTINATION lib )
install( FILES ${ROSALIND_LIB_HEADER_FILES} DESTINATION include/${ROSALIND_LIB} )


# Includes Catch in the project:
add_subdirectory( catch )
add_subdirectory( app )
add_subdirectory( test )




