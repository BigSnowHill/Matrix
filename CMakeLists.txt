cmake_minimum_required(VERSION 3.16)
project(lab_13)

set(CMAKE_CXX_STANDARD 14)

set(Header parallel_matrix_async.h)
set(Sources parallel_matrix_async.cpp)

add_compile_options(-pthread)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

add_library(parallel_matrix.h STATIC ${Sources} ${Headers})

add_executable(lab_13 main.cpp parallel_matrix_async.cpp parallel_matrix_async.h)



