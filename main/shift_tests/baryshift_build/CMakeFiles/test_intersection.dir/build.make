# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/twebb8/projects/cgal/c_curved_space/main/shift_tests/baryshift_source

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/twebb8/projects/cgal/c_curved_space/main/shift_tests/baryshift_build

# Include any dependencies generated for this target.
include CMakeFiles/test_intersection.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/test_intersection.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test_intersection.dir/flags.make

CMakeFiles/test_intersection.dir/baryshift.cpp.o: CMakeFiles/test_intersection.dir/flags.make
CMakeFiles/test_intersection.dir/baryshift.cpp.o: /home/twebb8/projects/cgal/c_curved_space/main/shift_tests/baryshift_source/baryshift.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/twebb8/projects/cgal/c_curved_space/main/shift_tests/baryshift_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/test_intersection.dir/baryshift.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_intersection.dir/baryshift.cpp.o -c /home/twebb8/projects/cgal/c_curved_space/main/shift_tests/baryshift_source/baryshift.cpp

CMakeFiles/test_intersection.dir/baryshift.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_intersection.dir/baryshift.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/twebb8/projects/cgal/c_curved_space/main/shift_tests/baryshift_source/baryshift.cpp > CMakeFiles/test_intersection.dir/baryshift.cpp.i

CMakeFiles/test_intersection.dir/baryshift.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_intersection.dir/baryshift.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/twebb8/projects/cgal/c_curved_space/main/shift_tests/baryshift_source/baryshift.cpp -o CMakeFiles/test_intersection.dir/baryshift.cpp.s

# Object files for target test_intersection
test_intersection_OBJECTS = \
"CMakeFiles/test_intersection.dir/baryshift.cpp.o"

# External object files for target test_intersection
test_intersection_EXTERNAL_OBJECTS =

test_intersection: CMakeFiles/test_intersection.dir/baryshift.cpp.o
test_intersection: CMakeFiles/test_intersection.dir/build.make
test_intersection: /usr/lib/x86_64-linux-gnu/libgmpxx.so
test_intersection: /usr/lib/x86_64-linux-gnu/libmpfr.so
test_intersection: /usr/lib/x86_64-linux-gnu/libgmp.so
test_intersection: CMakeFiles/test_intersection.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/twebb8/projects/cgal/c_curved_space/main/shift_tests/baryshift_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_intersection"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_intersection.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test_intersection.dir/build: test_intersection

.PHONY : CMakeFiles/test_intersection.dir/build

CMakeFiles/test_intersection.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test_intersection.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test_intersection.dir/clean

CMakeFiles/test_intersection.dir/depend:
	cd /home/twebb8/projects/cgal/c_curved_space/main/shift_tests/baryshift_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/twebb8/projects/cgal/c_curved_space/main/shift_tests/baryshift_source /home/twebb8/projects/cgal/c_curved_space/main/shift_tests/baryshift_source /home/twebb8/projects/cgal/c_curved_space/main/shift_tests/baryshift_build /home/twebb8/projects/cgal/c_curved_space/main/shift_tests/baryshift_build /home/twebb8/projects/cgal/c_curved_space/main/shift_tests/baryshift_build/CMakeFiles/test_intersection.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test_intersection.dir/depend

