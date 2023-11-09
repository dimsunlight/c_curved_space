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

add_executable( simulate  get_force.cpp shift.cpp utils.cpp simulate.cpp )
add_executable( force_test force_dist_test.cpp get_force.cpp ) 
add_executable( random_simulate get_force.cpp shift.cpp utils.cpp random_pos_simulate.cpp)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS simulate )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS force_test )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS random_simulate )

# Link the executable to CGAL and third-party libraries
target_link_libraries(simulate PRIVATE CGAL::CGAL )
target_link_libraries(force_test PRIVATE CGAL::CGAL )
target_link_libraries(random_simulate PRIVATE CGAL::CGAL ) 
