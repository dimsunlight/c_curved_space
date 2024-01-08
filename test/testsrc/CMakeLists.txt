cmake_minimum_required(VERSION 3.1...3.23)

add_library(src STATIC get_force.cpp shift.cpp)

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

target_include_directories(src PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
