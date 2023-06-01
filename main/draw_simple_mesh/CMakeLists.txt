cmake_minimum_required(VERSION 3.15)
project(draw_simple_mesh)

message("Generated with config types: ${CMAKE_CONFIGURATION_TYPES}")

set(CMAKE_BUILD_TYPE Release)

find_package(CGAL)
find_package(CGAL COMPONENTS Qt5)

if(CGAL_Qt5_FOUND)
    add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
endif()


add_executable(draw_simple_mesh draw_simple_mesh.cpp)

target_link_libraries(draw_simple_mesh PUBLIC  CGAL::CGAL)

if(CGAL_Qt5_FOUND)
  target_link_libraries(draw_triangulation_3 PUBLIC CGAL::CGAL_Qt5)
endif()
