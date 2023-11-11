cmake_minimum_required(VERSION 3.1...3.23)

project (shift) 

find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()

endif()

find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

add_executable(shift shift.cpp utils.cpp)
add_executable(test_intersection intersection_checker.cpp)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS shift)
add_to_cached_list( CGAL_EXECUTABLE_TARGETS test_intersection)

target_link_libraries(shift PRIVATE CGAL::CGAL)
target_link_libraries(test_intersection PRIVATE CGAL::CGAL)
