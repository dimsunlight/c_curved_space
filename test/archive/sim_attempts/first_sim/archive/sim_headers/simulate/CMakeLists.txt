# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project( simulate )


# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

create_single_source_cgal_program( "simulate.cpp" )

# include for local directory
add_subdirectory(force_and_shift)

target_link_libraries(simulate PUBLIC get_force)
target_link_libraries(simulate PUBLIC shift)

target_include_directories(simulate PUBLIC get_force)
target_include_directories(simulate PUBLIC shift)
# include for local package


# Creating entries for all C++ files with "main" routine
# ##########################################################






