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
include ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/depend.make

# Include the progress variables for this target.
include ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/progress.make

# Include the compile flags for this target's objects.
include ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/flags.make

ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/src/bitwise.c.o: ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/flags.make
ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/src/bitwise.c.o: ../ThirdParty/ogg/vtkogg/src/bitwise.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yyl/SEMO/lib/VTK/VTK-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/src/bitwise.c.o"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/ogg/vtkogg && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/ogg.dir/src/bitwise.c.o   -c /home/yyl/SEMO/lib/VTK/ThirdParty/ogg/vtkogg/src/bitwise.c

ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/src/bitwise.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/ogg.dir/src/bitwise.c.i"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/ogg/vtkogg && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/yyl/SEMO/lib/VTK/ThirdParty/ogg/vtkogg/src/bitwise.c > CMakeFiles/ogg.dir/src/bitwise.c.i

ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/src/bitwise.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/ogg.dir/src/bitwise.c.s"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/ogg/vtkogg && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/yyl/SEMO/lib/VTK/ThirdParty/ogg/vtkogg/src/bitwise.c -o CMakeFiles/ogg.dir/src/bitwise.c.s

ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/src/framing.c.o: ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/flags.make
ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/src/framing.c.o: ../ThirdParty/ogg/vtkogg/src/framing.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yyl/SEMO/lib/VTK/VTK-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/src/framing.c.o"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/ogg/vtkogg && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/ogg.dir/src/framing.c.o   -c /home/yyl/SEMO/lib/VTK/ThirdParty/ogg/vtkogg/src/framing.c

ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/src/framing.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/ogg.dir/src/framing.c.i"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/ogg/vtkogg && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/yyl/SEMO/lib/VTK/ThirdParty/ogg/vtkogg/src/framing.c > CMakeFiles/ogg.dir/src/framing.c.i

ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/src/framing.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/ogg.dir/src/framing.c.s"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/ogg/vtkogg && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/yyl/SEMO/lib/VTK/ThirdParty/ogg/vtkogg/src/framing.c -o CMakeFiles/ogg.dir/src/framing.c.s

# Object files for target ogg
ogg_OBJECTS = \
"CMakeFiles/ogg.dir/src/bitwise.c.o" \
"CMakeFiles/ogg.dir/src/framing.c.o"

# External object files for target ogg
ogg_EXTERNAL_OBJECTS =

lib/libvtkogg-9.0.so.9.0.1: ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/src/bitwise.c.o
lib/libvtkogg-9.0.so.9.0.1: ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/src/framing.c.o
lib/libvtkogg-9.0.so.9.0.1: ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/build.make
lib/libvtkogg-9.0.so.9.0.1: ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/yyl/SEMO/lib/VTK/VTK-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C shared library ../../../lib/libvtkogg-9.0.so"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/ogg/vtkogg && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ogg.dir/link.txt --verbose=$(VERBOSE)
	cd /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/ogg/vtkogg && $(CMAKE_COMMAND) -E cmake_symlink_library ../../../lib/libvtkogg-9.0.so.9.0.1 ../../../lib/libvtkogg-9.0.so.1 ../../../lib/libvtkogg-9.0.so

lib/libvtkogg-9.0.so.1: lib/libvtkogg-9.0.so.9.0.1
	@$(CMAKE_COMMAND) -E touch_nocreate lib/libvtkogg-9.0.so.1

lib/libvtkogg-9.0.so: lib/libvtkogg-9.0.so.9.0.1
	@$(CMAKE_COMMAND) -E touch_nocreate lib/libvtkogg-9.0.so

# Rule to build all files generated by this target.
ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/build: lib/libvtkogg-9.0.so

.PHONY : ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/build

ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/clean:
	cd /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/ogg/vtkogg && $(CMAKE_COMMAND) -P CMakeFiles/ogg.dir/cmake_clean.cmake
.PHONY : ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/clean

ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/depend:
	cd /home/yyl/SEMO/lib/VTK/VTK-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yyl/SEMO/lib/VTK /home/yyl/SEMO/lib/VTK/ThirdParty/ogg/vtkogg /home/yyl/SEMO/lib/VTK/VTK-build /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/ogg/vtkogg /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ThirdParty/ogg/vtkogg/CMakeFiles/ogg.dir/depend

