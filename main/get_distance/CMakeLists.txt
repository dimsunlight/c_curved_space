# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Surface_mesh_shortest_path_Examples )



find_package(CGAL QUIET)

if ( CGAL_FOUND )

  create_single_source_cgal_program( "multipoint_distance.cpp" )
  create_single_source_cgal_program( "shortest_paths.cpp")
  
else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()


