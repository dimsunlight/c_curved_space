cmake_minimum_required(VERSION 3.23)

project(Simple_Mesh)

find_package(CGAL)
add_executable(simple_mesh simple_mesh.cpp)
target_link_libraries(simple_mesh CGAL::CGAL)
