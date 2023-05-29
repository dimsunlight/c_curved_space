cmake_minimum_required(VERSION 3.20)

project(Simple_Mesh)

message("Generated with config types: ${CMAKE_CONFIGURATION_TYPES}")

set(CMAKE_BUILD_TYPE Release)

find_package(CGAL)
add_executable(simple_mesh simple_mesh.cpp)
target_link_libraries(simple_mesh CGAL::CGAL)
