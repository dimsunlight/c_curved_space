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
CMAKE_SOURCE_DIR = /home/twebb8/projects/cgal/c_curved_space/main/attempt_simulate/sim_new/sim_new_source

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/twebb8/projects/cgal/c_curved_space/main/attempt_simulate/sim_new/sim_new_build

# Include any dependencies generated for this target.
include CMakeFiles/pathbend_shift.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/pathbend_shift.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/pathbend_shift.dir/flags.make

CMakeFiles/pathbend_shift.dir/pathbend_shift.cpp.o: CMakeFiles/pathbend_shift.dir/flags.make
CMakeFiles/pathbend_shift.dir/pathbend_shift.cpp.o: /home/twebb8/projects/cgal/c_curved_space/main/attempt_simulate/sim_new/sim_new_source/pathbend_shift.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/twebb8/projects/cgal/c_curved_space/main/attempt_simulate/sim_new/sim_new_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/pathbend_shift.dir/pathbend_shift.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/pathbend_shift.dir/pathbend_shift.cpp.o -c /home/twebb8/projects/cgal/c_curved_space/main/attempt_simulate/sim_new/sim_new_source/pathbend_shift.cpp

CMakeFiles/pathbend_shift.dir/pathbend_shift.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pathbend_shift.dir/pathbend_shift.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/twebb8/projects/cgal/c_curved_space/main/attempt_simulate/sim_new/sim_new_source/pathbend_shift.cpp > CMakeFiles/pathbend_shift.dir/pathbend_shift.cpp.i

CMakeFiles/pathbend_shift.dir/pathbend_shift.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pathbend_shift.dir/pathbend_shift.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/twebb8/projects/cgal/c_curved_space/main/attempt_simulate/sim_new/sim_new_source/pathbend_shift.cpp -o CMakeFiles/pathbend_shift.dir/pathbend_shift.cpp.s

# Object files for target pathbend_shift
pathbend_shift_OBJECTS = \
"CMakeFiles/pathbend_shift.dir/pathbend_shift.cpp.o"

# External object files for target pathbend_shift
pathbend_shift_EXTERNAL_OBJECTS =

pathbend_shift: CMakeFiles/pathbend_shift.dir/pathbend_shift.cpp.o
pathbend_shift: CMakeFiles/pathbend_shift.dir/build.make
pathbend_shift: /usr/lib/x86_64-linux-gnu/libgmpxx.so
pathbend_shift: /usr/lib/x86_64-linux-gnu/libmpfr.so
pathbend_shift: /usr/lib/x86_64-linux-gnu/libgmp.so
pathbend_shift: CMakeFiles/pathbend_shift.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/twebb8/projects/cgal/c_curved_space/main/attempt_simulate/sim_new/sim_new_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable pathbend_shift"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pathbend_shift.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/pathbend_shift.dir/build: pathbend_shift

.PHONY : CMakeFiles/pathbend_shift.dir/build

CMakeFiles/pathbend_shift.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/pathbend_shift.dir/cmake_clean.cmake
.PHONY : CMakeFiles/pathbend_shift.dir/clean

CMakeFiles/pathbend_shift.dir/depend:
	cd /home/twebb8/projects/cgal/c_curved_space/main/attempt_simulate/sim_new/sim_new_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/twebb8/projects/cgal/c_curved_space/main/attempt_simulate/sim_new/sim_new_source /home/twebb8/projects/cgal/c_curved_space/main/attempt_simulate/sim_new/sim_new_source /home/twebb8/projects/cgal/c_curved_space/main/attempt_simulate/sim_new/sim_new_build /home/twebb8/projects/cgal/c_curved_space/main/attempt_simulate/sim_new/sim_new_build /home/twebb8/projects/cgal/c_curved_space/main/attempt_simulate/sim_new/sim_new_build/CMakeFiles/pathbend_shift.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/pathbend_shift.dir/depend

