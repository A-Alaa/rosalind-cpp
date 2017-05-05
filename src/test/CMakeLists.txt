# Make test executable


configure_file( tests.h.in ${CMAKE_CURRENT_SOURCE_DIR}/tests.h @ONLY)


set(TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests.cpp)
set(TEST_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/tests.h)
add_executable(tests ${TEST_SOURCES} ${TEST_HEADERS})
target_link_libraries(tests Catch)
