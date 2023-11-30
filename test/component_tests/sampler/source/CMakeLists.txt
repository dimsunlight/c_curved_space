cmake_minimum_required(VERSION 3.1...3.23)

project( force )

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

add_executable( sp_sample sp_sample.cpp utils.cpp ) 
add_executable( torus_sample tor_sample.cpp utils.cpp ) 

add_to_cached_list( CGAL_EXECUTABLE_TARGETS sp_sample )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS torus_sample )

target_link_libraries(sp_sample PRIVATE CGAL::CGAL ) 
target_link_libraries(torus_sample PRIVATE CGAL::CGAL )  
