cmake_minimum_required(VERSION 3.13)

project(base)


add_subdirectory(function)
add_executable(base base.cpp)
target_link_libraries(base PUBLIC function)
target_include_directories(base PUBLIC function)
