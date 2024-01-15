cmake_minimum_required(VERSION 3.1...3.28)
project(
  gprof_test
  VERSION 1.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pg") #gprof profiling flags

add_executable(gprof_hw gprof_hw.cpp)

