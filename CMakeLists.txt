cmake_minimum_required(VERSION 3.10)
project(GANGSTA_plus)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags equivalent to the existing Makefile
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -fprefetch-loop-arrays -mfpmath=sse -funswitch-loops -DNDEBUG -Wno-deprecated")

# Add the executable
add_executable(gplus gplus.cpp)

# Install target
install(TARGETS gplus DESTINATION bin)

# Output message for clarity
message(STATUS "Build with: cmake -B build && cmake --build build")
message(STATUS "Run with: ./build/gplus d2uaga1.pdb d1gkub1.pdb")