cmake_minimum_required(VERSION 3.0)

project(CTRU-Compact)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c99")

add_compile_options(
    -Wall 
    -march=native 
    -mtune=native 
    -O3 
    -fomit-frame-pointer 
    -Wno-unknown-pragmas
    -lcrypto 
    -ldl)

add_subdirectory(common)
add_subdirectory(src)
add_subdirectory(test)
