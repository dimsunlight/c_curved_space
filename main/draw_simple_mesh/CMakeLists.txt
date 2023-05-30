cmake_minimum_required(VERSION 3.20)
project(draw_simple_mesh)

message("Generated with config types: ${CMAKE_CONFIGURATION_TYPES}")

set(CMAKE_BUILD_TYPE Release)

find_package(CGAL)
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

target_link_libraries(draw_simple_mesh CGAL::CGAL)

if(CGAL_Qt5_FOUND)
  #link it with the required CGAL libraries
  target_link_libraries(draw_simple_mesh PUBLIC CGAL::CGAL_Basic_viewer)
  target_link_libraries(sm_draw_small_faces PUBLIC CGAL::CGAL_Basic_viewer)
endif()

add_executable(draw_simple_mesh draw_simple_mesh.cpp)
