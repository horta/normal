cmake_minimum_required (VERSION 3.4)
project (normal)
set(CMAKE_BUILD_TYPE Release)

include_directories(.)
add_subdirectory (test)

enable_testing()
add_test( NAME TestFunctions COMMAND test_functions )
