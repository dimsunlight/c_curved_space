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
CMAKE_SOURCE_DIR = /home/toler/projects/CGAL_curvature/main/get_distance

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/toler/projects/CGAL_curvature/main/get_distance_build

# Include any dependencies generated for this target.
include CMakeFiles/using_locate.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/using_locate.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/using_locate.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/using_locate.dir/flags.make

CMakeFiles/using_locate.dir/using_locate.cpp.o: CMakeFiles/using_locate.dir/flags.make
CMakeFiles/using_locate.dir/using_locate.cpp.o: /home/toler/projects/CGAL_curvature/main/get_distance/using_locate.cpp
CMakeFiles/using_locate.dir/using_locate.cpp.o: CMakeFiles/using_locate.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/toler/projects/CGAL_curvature/main/get_distance_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/using_locate.dir/using_locate.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/using_locate.dir/using_locate.cpp.o -MF CMakeFiles/using_locate.dir/using_locate.cpp.o.d -o CMakeFiles/using_locate.dir/using_locate.cpp.o -c /home/toler/projects/CGAL_curvature/main/get_distance/using_locate.cpp

CMakeFiles/using_locate.dir/using_locate.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/using_locate.dir/using_locate.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/toler/projects/CGAL_curvature/main/get_distance/using_locate.cpp > CMakeFiles/using_locate.dir/using_locate.cpp.i

CMakeFiles/using_locate.dir/using_locate.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/using_locate.dir/using_locate.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/toler/projects/CGAL_curvature/main/get_distance/using_locate.cpp -o CMakeFiles/using_locate.dir/using_locate.cpp.s

# Object files for target using_locate
using_locate_OBJECTS = \
"CMakeFiles/using_locate.dir/using_locate.cpp.o"

# External object files for target using_locate
using_locate_EXTERNAL_OBJECTS =

using_locate: CMakeFiles/using_locate.dir/using_locate.cpp.o
using_locate: CMakeFiles/using_locate.dir/build.make
using_locate: /usr/lib/x86_64-linux-gnu/libgmpxx.so
using_locate: /usr/lib/x86_64-linux-gnu/libmpfr.so
using_locate: /usr/lib/x86_64-linux-gnu/libgmp.so
using_locate: CMakeFiles/using_locate.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/toler/projects/CGAL_curvature/main/get_distance_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable using_locate"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/using_locate.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/using_locate.dir/build: using_locate
.PHONY : CMakeFiles/using_locate.dir/build

CMakeFiles/using_locate.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/using_locate.dir/cmake_clean.cmake
.PHONY : CMakeFiles/using_locate.dir/clean

CMakeFiles/using_locate.dir/depend:
	cd /home/toler/projects/CGAL_curvature/main/get_distance_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/toler/projects/CGAL_curvature/main/get_distance /home/toler/projects/CGAL_curvature/main/get_distance /home/toler/projects/CGAL_curvature/main/get_distance_build /home/toler/projects/CGAL_curvature/main/get_distance_build /home/toler/projects/CGAL_curvature/main/get_distance_build/CMakeFiles/using_locate.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/using_locate.dir/depend

