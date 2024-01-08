cmake_minimum_required(VERSION 3.0.0)

project(ccs_test LANGUAGES cxx)

add_subdirectory(testsrc)
add_subdirectory(utils)

add_subdirectory(component_tests)
add_subdirectory(mesher) 
add_subdirectory(simulate) 

