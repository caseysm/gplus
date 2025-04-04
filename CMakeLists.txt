cmake_minimum_required(VERSION 3.10)
project(GANGSTA_plus VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags equivalent to the existing Makefile
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -fprefetch-loop-arrays -mfpmath=sse -funswitch-loops -DNDEBUG -Wno-deprecated")

# Define directories
set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(APPS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/apps")
set(EXAMPLES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/examples")
set(TESTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests")

# Set up include directories
include_directories(${INCLUDE_DIR})

# Collect all source files
file(GLOB_RECURSE SOURCES 
    "${SRC_DIR}/*.cpp"
)

# Library target
add_library(gangsta_lib STATIC ${SOURCES})

# Main executable
add_executable(gplus "${APPS_DIR}/gplus.cpp")
target_link_libraries(gplus PRIVATE gangsta_lib)

# Copy example PDB files to build directory
file(GLOB PDB_FILES "${EXAMPLES_DIR}/*.pdb")
foreach(PDB_FILE ${PDB_FILES})
    file(COPY ${PDB_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

# Install targets
install(TARGETS gplus DESTINATION bin)

# No external dependencies required

# Testing with Google Test
enable_testing()

# Download and configure Google Test
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/v1.13.0.tar.gz
)
FetchContent_MakeAvailable(googletest)

# Include Google Test headers
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Scan for test files (excluding debug tool)
file(GLOB_RECURSE TEST_SOURCES "${TESTS_DIR}/*.cpp")
list(REMOVE_ITEM TEST_SOURCES "${TESTS_DIR}/algorithm/debug_alignment.cpp")
if(TEST_SOURCES)
    add_executable(run_tests ${TEST_SOURCES})
    target_link_libraries(run_tests PRIVATE 
        gangsta_lib 
        gtest 
        gtest_main
        pthread
    )
    add_test(NAME AllTests COMMAND run_tests)
endif()

# Helper function to add individual test executables
function(add_unit_test test_name)
    add_executable(${test_name} ${TESTS_DIR}/${test_name}.cpp)
    target_link_libraries(${test_name} PRIVATE 
        gangsta_lib 
        gtest 
        gtest_main
        pthread
    )
    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

# Add debug alignment tool
add_executable(debug_alignment "${TESTS_DIR}/algorithm/debug_alignment.cpp")
target_link_libraries(debug_alignment PRIVATE gangsta_lib)

# Output message for clarity
message(STATUS "Build with: cmake -B build && cmake --build build")
message(STATUS "Run with: ./build/gplus d2uaga1.pdb d1gkub1.pdb")