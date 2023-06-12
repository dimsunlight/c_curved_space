# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project( simple_mesh )

# CGAL and its components

find_package(CGAL)

create_single_source_cgal_program("heightmap_mesh.cpp")
create_single_source_cgal_program("torus_mesh.cpp")
create_single_source_cgal_program("compute_normals.cpp")
