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
CMAKE_SOURCE_DIR = /home/yyl/SEMO/lib/VTK

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/yyl/SEMO/lib/VTK/VTK-build

# Include any dependencies generated for this target.
include Wrapping/Tools/CMakeFiles/WrapJava.dir/depend.make

# Include the progress variables for this target.
include Wrapping/Tools/CMakeFiles/WrapJava.dir/progress.make

# Include the compile flags for this target's objects.
include Wrapping/Tools/CMakeFiles/WrapJava.dir/flags.make

Wrapping/Tools/CMakeFiles/WrapJava.dir/vtkWrapJava.c.o: Wrapping/Tools/CMakeFiles/WrapJava.dir/flags.make
Wrapping/Tools/CMakeFiles/WrapJava.dir/vtkWrapJava.c.o: ../Wrapping/Tools/vtkWrapJava.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yyl/SEMO/lib/VTK/VTK-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object Wrapping/Tools/CMakeFiles/WrapJava.dir/vtkWrapJava.c.o"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/Wrapping/Tools && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/WrapJava.dir/vtkWrapJava.c.o   -c /home/yyl/SEMO/lib/VTK/Wrapping/Tools/vtkWrapJava.c

Wrapping/Tools/CMakeFiles/WrapJava.dir/vtkWrapJava.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/WrapJava.dir/vtkWrapJava.c.i"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/Wrapping/Tools && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/yyl/SEMO/lib/VTK/Wrapping/Tools/vtkWrapJava.c > CMakeFiles/WrapJava.dir/vtkWrapJava.c.i

Wrapping/Tools/CMakeFiles/WrapJava.dir/vtkWrapJava.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/WrapJava.dir/vtkWrapJava.c.s"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/Wrapping/Tools && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/yyl/SEMO/lib/VTK/Wrapping/Tools/vtkWrapJava.c -o CMakeFiles/WrapJava.dir/vtkWrapJava.c.s

# Object files for target WrapJava
WrapJava_OBJECTS = \
"CMakeFiles/WrapJava.dir/vtkWrapJava.c.o"

# External object files for target WrapJava
WrapJava_EXTERNAL_OBJECTS =

bin/vtkWrapJava-9.0: Wrapping/Tools/CMakeFiles/WrapJava.dir/vtkWrapJava.c.o
bin/vtkWrapJava-9.0: Wrapping/Tools/CMakeFiles/WrapJava.dir/build.make
bin/vtkWrapJava-9.0: lib/libvtkWrappingTools-9.0.so.9.0.1
bin/vtkWrapJava-9.0: Wrapping/Tools/CMakeFiles/WrapJava.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/yyl/SEMO/lib/VTK/VTK-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable ../../bin/vtkWrapJava-9.0"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/Wrapping/Tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/WrapJava.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Wrapping/Tools/CMakeFiles/WrapJava.dir/build: bin/vtkWrapJava-9.0

.PHONY : Wrapping/Tools/CMakeFiles/WrapJava.dir/build

Wrapping/Tools/CMakeFiles/WrapJava.dir/clean:
	cd /home/yyl/SEMO/lib/VTK/VTK-build/Wrapping/Tools && $(CMAKE_COMMAND) -P CMakeFiles/WrapJava.dir/cmake_clean.cmake
.PHONY : Wrapping/Tools/CMakeFiles/WrapJava.dir/clean

Wrapping/Tools/CMakeFiles/WrapJava.dir/depend:
	cd /home/yyl/SEMO/lib/VTK/VTK-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yyl/SEMO/lib/VTK /home/yyl/SEMO/lib/VTK/Wrapping/Tools /home/yyl/SEMO/lib/VTK/VTK-build /home/yyl/SEMO/lib/VTK/VTK-build/Wrapping/Tools /home/yyl/SEMO/lib/VTK/VTK-build/Wrapping/Tools/CMakeFiles/WrapJava.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Wrapping/Tools/CMakeFiles/WrapJava.dir/depend

