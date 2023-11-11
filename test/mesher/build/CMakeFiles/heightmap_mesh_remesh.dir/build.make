# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/toler/projects/CGAL_curvature/main/simple_mesh

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/toler/projects/CGAL_curvature/main/simple_mesh_build

# Include any dependencies generated for this target.
include CMakeFiles/heightmap_mesh_remesh.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/heightmap_mesh_remesh.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/heightmap_mesh_remesh.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/heightmap_mesh_remesh.dir/flags.make

CMakeFiles/heightmap_mesh_remesh.dir/heightmap_mesh_remesh.cpp.o: CMakeFiles/heightmap_mesh_remesh.dir/flags.make
CMakeFiles/heightmap_mesh_remesh.dir/heightmap_mesh_remesh.cpp.o: /home/toler/projects/CGAL_curvature/main/simple_mesh/heightmap_mesh_remesh.cpp
CMakeFiles/heightmap_mesh_remesh.dir/heightmap_mesh_remesh.cpp.o: CMakeFiles/heightmap_mesh_remesh.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/toler/projects/CGAL_curvature/main/simple_mesh_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/heightmap_mesh_remesh.dir/heightmap_mesh_remesh.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/heightmap_mesh_remesh.dir/heightmap_mesh_remesh.cpp.o -MF CMakeFiles/heightmap_mesh_remesh.dir/heightmap_mesh_remesh.cpp.o.d -o CMakeFiles/heightmap_mesh_remesh.dir/heightmap_mesh_remesh.cpp.o -c /home/toler/projects/CGAL_curvature/main/simple_mesh/heightmap_mesh_remesh.cpp

CMakeFiles/heightmap_mesh_remesh.dir/heightmap_mesh_remesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/heightmap_mesh_remesh.dir/heightmap_mesh_remesh.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/toler/projects/CGAL_curvature/main/simple_mesh/heightmap_mesh_remesh.cpp > CMakeFiles/heightmap_mesh_remesh.dir/heightmap_mesh_remesh.cpp.i

CMakeFiles/heightmap_mesh_remesh.dir/heightmap_mesh_remesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/heightmap_mesh_remesh.dir/heightmap_mesh_remesh.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/toler/projects/CGAL_curvature/main/simple_mesh/heightmap_mesh_remesh.cpp -o CMakeFiles/heightmap_mesh_remesh.dir/heightmap_mesh_remesh.cpp.s

# Object files for target heightmap_mesh_remesh
heightmap_mesh_remesh_OBJECTS = \
"CMakeFiles/heightmap_mesh_remesh.dir/heightmap_mesh_remesh.cpp.o"

# External object files for target heightmap_mesh_remesh
heightmap_mesh_remesh_EXTERNAL_OBJECTS =

heightmap_mesh_remesh: CMakeFiles/heightmap_mesh_remesh.dir/heightmap_mesh_remesh.cpp.o
heightmap_mesh_remesh: CMakeFiles/heightmap_mesh_remesh.dir/build.make
heightmap_mesh_remesh: /usr/lib/x86_64-linux-gnu/libgmpxx.so
heightmap_mesh_remesh: /usr/lib/x86_64-linux-gnu/libmpfr.so
heightmap_mesh_remesh: /usr/lib/x86_64-linux-gnu/libgmp.so
heightmap_mesh_remesh: CMakeFiles/heightmap_mesh_remesh.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/toler/projects/CGAL_curvature/main/simple_mesh_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable heightmap_mesh_remesh"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/heightmap_mesh_remesh.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/heightmap_mesh_remesh.dir/build: heightmap_mesh_remesh
.PHONY : CMakeFiles/heightmap_mesh_remesh.dir/build

CMakeFiles/heightmap_mesh_remesh.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/heightmap_mesh_remesh.dir/cmake_clean.cmake
.PHONY : CMakeFiles/heightmap_mesh_remesh.dir/clean

CMakeFiles/heightmap_mesh_remesh.dir/depend:
	cd /home/toler/projects/CGAL_curvature/main/simple_mesh_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/toler/projects/CGAL_curvature/main/simple_mesh /home/toler/projects/CGAL_curvature/main/simple_mesh /home/toler/projects/CGAL_curvature/main/simple_mesh_build /home/toler/projects/CGAL_curvature/main/simple_mesh_build /home/toler/projects/CGAL_curvature/main/simple_mesh_build/CMakeFiles/heightmap_mesh_remesh.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/heightmap_mesh_remesh.dir/depend

