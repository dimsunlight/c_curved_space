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

# include for local directory

# include for local package


# Creating entries for target: simulate
# ############################

add_executable( force_test force_dist_test.cpp get_force.cpp utils.cpp ) 
add_executable( f_v_n_test force_per_target_test.cpp get_force.cpp utils.cpp ) 
add_executable( seq_test seq_test.cpp get_force.cpp utils.cpp ) 

add_to_cached_list( CGAL_EXECUTABLE_TARGETS force_test )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS f_v_n_test ) 
add_to_cached_list( CGAL_EXECUTABLE_TARGETS seq_test ) 

# Link the executable to CGAL and third-party libraries
target_link_libraries(force_test PRIVATE CGAL::CGAL )
target_link_libraries(f_v_n_test PRIVATE CGAL::CGAL ) 
target_link_libraries(seq_test PRIVATE CGAL::CGAL ) 
