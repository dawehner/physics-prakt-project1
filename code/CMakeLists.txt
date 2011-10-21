project(nbody)
cmake_minimum_required(VERSION 2.8)
# set(CMAKE_CXX_COMPILER "clang++")
# set(CMAKE_LINKER "llvm-ld-2.8")
set(CMAKE_CXX_FLAGS "-g -Wall -o3 -std=gnu++0x")

FIND_PACKAGE( Boost 1.40 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

add_executable(nbody main.cpp)

target_link_libraries(nbody ${Boost_PROGRAM_OPTIONS_LIBRARY} )
