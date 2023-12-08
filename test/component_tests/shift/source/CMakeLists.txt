# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project( shift )


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


# Creating entries for target: shift
# ############################

add_executable( shift_test shift_test.cpp shift.cpp utils.cpp )  
add_executable( vert_test vert_test.cpp shift.cpp utils.cpp ) 

add_to_cached_list( CGAL_EXECUTABLE_TARGETS shift_test ) 
add_to_cached_list( CGAL_EXECUTABLE_TARGETS vert_test ) 

# Link the executable to CGAL and third-party libraries

target_link_libraries(shift_test PRIVATE CGAL::CGAL )

target_link_libraries(vert_test PRIVATE CGAL::CGAL ) 

