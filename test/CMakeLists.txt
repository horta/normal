include_directories (..)
file(GLOB FILES "*.txt")
file(COPY ${FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_executable (test_functions test_functions.c)
