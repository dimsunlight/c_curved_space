# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project( baryshift_source )


# CGAL and its components

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  create_single_source_cgal_program("baryshift.cpp")

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

