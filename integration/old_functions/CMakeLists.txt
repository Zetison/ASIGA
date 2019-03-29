cmake_minimum_required(VERSION 2.6)
project(exactSol)

find_package(MPI)
enable_language(C)

include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_C_FLAGS "-std=c99")
add_executable(exactSol exactSol.c gaussQuad.h)

target_link_libraries(exactSol m)
