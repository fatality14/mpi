cmake_minimum_required(VERSION 3.5)

project(mpii LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(HEADER_FILES headers/mpi.h)

add_executable(mpii main.cpp ${HEADER_FILES})
target_compile_options(mpii PRIVATE -Wall -Wextra -pedantic)

find_library(MSMPI msmpi.lib libs)
target_link_libraries(mpii LINK_PUBLIC ${MSMPI})
