cmake_minimum_required(VERSION 3.11...3.23)

set(CMAKE_BUILD_TYPE "Release")

find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/testsrc
	${CMAKE_CURRENT_SOURCE_DIR}/util
   )

add_library(src STATIC get_force.cpp shift.cpp)
add_library(util STATIC utils.cpp)


add_subdirectory(testsrc)
add_subdirectory(utils)

add_subdirectory(component_tests)
add_subdirectory(mesher) 
add_subdirectory(simulate) 

