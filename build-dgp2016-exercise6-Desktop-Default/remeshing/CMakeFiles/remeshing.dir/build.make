# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_SOURCE_DIR = /home/jad/Documents/DigitalGeralt-/dgp2016-exercise6

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default

# Include any dependencies generated for this target.
include remeshing/CMakeFiles/remeshing.dir/depend.make

# Include the progress variables for this target.
include remeshing/CMakeFiles/remeshing.dir/progress.make

# Include the compile flags for this target's objects.
include remeshing/CMakeFiles/remeshing.dir/flags.make

remeshing/CMakeFiles/remeshing.dir/main.cpp.o: remeshing/CMakeFiles/remeshing.dir/flags.make
remeshing/CMakeFiles/remeshing.dir/main.cpp.o: /home/jad/Documents/DigitalGeralt-/dgp2016-exercise6/remeshing/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object remeshing/CMakeFiles/remeshing.dir/main.cpp.o"
	cd /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/remeshing && /usr/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/remeshing.dir/main.cpp.o -c /home/jad/Documents/DigitalGeralt-/dgp2016-exercise6/remeshing/main.cpp

remeshing/CMakeFiles/remeshing.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/remeshing.dir/main.cpp.i"
	cd /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/remeshing && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jad/Documents/DigitalGeralt-/dgp2016-exercise6/remeshing/main.cpp > CMakeFiles/remeshing.dir/main.cpp.i

remeshing/CMakeFiles/remeshing.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/remeshing.dir/main.cpp.s"
	cd /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/remeshing && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jad/Documents/DigitalGeralt-/dgp2016-exercise6/remeshing/main.cpp -o CMakeFiles/remeshing.dir/main.cpp.s

remeshing/CMakeFiles/remeshing.dir/main.cpp.o.requires:

.PHONY : remeshing/CMakeFiles/remeshing.dir/main.cpp.o.requires

remeshing/CMakeFiles/remeshing.dir/main.cpp.o.provides: remeshing/CMakeFiles/remeshing.dir/main.cpp.o.requires
	$(MAKE) -f remeshing/CMakeFiles/remeshing.dir/build.make remeshing/CMakeFiles/remeshing.dir/main.cpp.o.provides.build
.PHONY : remeshing/CMakeFiles/remeshing.dir/main.cpp.o.provides

remeshing/CMakeFiles/remeshing.dir/main.cpp.o.provides.build: remeshing/CMakeFiles/remeshing.dir/main.cpp.o


remeshing/CMakeFiles/remeshing.dir/viewer.cpp.o: remeshing/CMakeFiles/remeshing.dir/flags.make
remeshing/CMakeFiles/remeshing.dir/viewer.cpp.o: /home/jad/Documents/DigitalGeralt-/dgp2016-exercise6/remeshing/viewer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object remeshing/CMakeFiles/remeshing.dir/viewer.cpp.o"
	cd /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/remeshing && /usr/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/remeshing.dir/viewer.cpp.o -c /home/jad/Documents/DigitalGeralt-/dgp2016-exercise6/remeshing/viewer.cpp

remeshing/CMakeFiles/remeshing.dir/viewer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/remeshing.dir/viewer.cpp.i"
	cd /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/remeshing && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jad/Documents/DigitalGeralt-/dgp2016-exercise6/remeshing/viewer.cpp > CMakeFiles/remeshing.dir/viewer.cpp.i

remeshing/CMakeFiles/remeshing.dir/viewer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/remeshing.dir/viewer.cpp.s"
	cd /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/remeshing && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jad/Documents/DigitalGeralt-/dgp2016-exercise6/remeshing/viewer.cpp -o CMakeFiles/remeshing.dir/viewer.cpp.s

remeshing/CMakeFiles/remeshing.dir/viewer.cpp.o.requires:

.PHONY : remeshing/CMakeFiles/remeshing.dir/viewer.cpp.o.requires

remeshing/CMakeFiles/remeshing.dir/viewer.cpp.o.provides: remeshing/CMakeFiles/remeshing.dir/viewer.cpp.o.requires
	$(MAKE) -f remeshing/CMakeFiles/remeshing.dir/build.make remeshing/CMakeFiles/remeshing.dir/viewer.cpp.o.provides.build
.PHONY : remeshing/CMakeFiles/remeshing.dir/viewer.cpp.o.provides

remeshing/CMakeFiles/remeshing.dir/viewer.cpp.o.provides.build: remeshing/CMakeFiles/remeshing.dir/viewer.cpp.o


remeshing/CMakeFiles/remeshing.dir/mesh_processing.cpp.o: remeshing/CMakeFiles/remeshing.dir/flags.make
remeshing/CMakeFiles/remeshing.dir/mesh_processing.cpp.o: /home/jad/Documents/DigitalGeralt-/dgp2016-exercise6/remeshing/mesh_processing.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object remeshing/CMakeFiles/remeshing.dir/mesh_processing.cpp.o"
	cd /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/remeshing && /usr/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/remeshing.dir/mesh_processing.cpp.o -c /home/jad/Documents/DigitalGeralt-/dgp2016-exercise6/remeshing/mesh_processing.cpp

remeshing/CMakeFiles/remeshing.dir/mesh_processing.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/remeshing.dir/mesh_processing.cpp.i"
	cd /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/remeshing && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jad/Documents/DigitalGeralt-/dgp2016-exercise6/remeshing/mesh_processing.cpp > CMakeFiles/remeshing.dir/mesh_processing.cpp.i

remeshing/CMakeFiles/remeshing.dir/mesh_processing.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/remeshing.dir/mesh_processing.cpp.s"
	cd /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/remeshing && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jad/Documents/DigitalGeralt-/dgp2016-exercise6/remeshing/mesh_processing.cpp -o CMakeFiles/remeshing.dir/mesh_processing.cpp.s

remeshing/CMakeFiles/remeshing.dir/mesh_processing.cpp.o.requires:

.PHONY : remeshing/CMakeFiles/remeshing.dir/mesh_processing.cpp.o.requires

remeshing/CMakeFiles/remeshing.dir/mesh_processing.cpp.o.provides: remeshing/CMakeFiles/remeshing.dir/mesh_processing.cpp.o.requires
	$(MAKE) -f remeshing/CMakeFiles/remeshing.dir/build.make remeshing/CMakeFiles/remeshing.dir/mesh_processing.cpp.o.provides.build
.PHONY : remeshing/CMakeFiles/remeshing.dir/mesh_processing.cpp.o.provides

remeshing/CMakeFiles/remeshing.dir/mesh_processing.cpp.o.provides.build: remeshing/CMakeFiles/remeshing.dir/mesh_processing.cpp.o


# Object files for target remeshing
remeshing_OBJECTS = \
"CMakeFiles/remeshing.dir/main.cpp.o" \
"CMakeFiles/remeshing.dir/viewer.cpp.o" \
"CMakeFiles/remeshing.dir/mesh_processing.cpp.o"

# External object files for target remeshing
remeshing_EXTERNAL_OBJECTS =

remeshing/remeshing: remeshing/CMakeFiles/remeshing.dir/main.cpp.o
remeshing/remeshing: remeshing/CMakeFiles/remeshing.dir/viewer.cpp.o
remeshing/remeshing: remeshing/CMakeFiles/remeshing.dir/mesh_processing.cpp.o
remeshing/remeshing: remeshing/CMakeFiles/remeshing.dir/build.make
remeshing/remeshing: externals/surface_mesh/libsurface_mesh.so.1.0
remeshing/remeshing: nanogui/libnanogui.a
remeshing/remeshing: remeshing/CMakeFiles/remeshing.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable remeshing"
	cd /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/remeshing && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/remeshing.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
remeshing/CMakeFiles/remeshing.dir/build: remeshing/remeshing

.PHONY : remeshing/CMakeFiles/remeshing.dir/build

remeshing/CMakeFiles/remeshing.dir/requires: remeshing/CMakeFiles/remeshing.dir/main.cpp.o.requires
remeshing/CMakeFiles/remeshing.dir/requires: remeshing/CMakeFiles/remeshing.dir/viewer.cpp.o.requires
remeshing/CMakeFiles/remeshing.dir/requires: remeshing/CMakeFiles/remeshing.dir/mesh_processing.cpp.o.requires

.PHONY : remeshing/CMakeFiles/remeshing.dir/requires

remeshing/CMakeFiles/remeshing.dir/clean:
	cd /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/remeshing && $(CMAKE_COMMAND) -P CMakeFiles/remeshing.dir/cmake_clean.cmake
.PHONY : remeshing/CMakeFiles/remeshing.dir/clean

remeshing/CMakeFiles/remeshing.dir/depend:
	cd /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jad/Documents/DigitalGeralt-/dgp2016-exercise6 /home/jad/Documents/DigitalGeralt-/dgp2016-exercise6/remeshing /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/remeshing /home/jad/Documents/DigitalGeralt-/build-dgp2016-exercise6-Desktop-Default/remeshing/CMakeFiles/remeshing.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : remeshing/CMakeFiles/remeshing.dir/depend
