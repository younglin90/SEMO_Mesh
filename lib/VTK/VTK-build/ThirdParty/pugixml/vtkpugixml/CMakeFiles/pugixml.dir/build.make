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
include ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/depend.make

# Include the progress variables for this target.
include ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/progress.make

# Include the compile flags for this target's objects.
include ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/flags.make

ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/src/pugixml.cpp.o: ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/flags.make
ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/src/pugixml.cpp.o: ../ThirdParty/pugixml/vtkpugixml/src/pugixml.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yyl/SEMO/lib/VTK/VTK-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/src/pugixml.cpp.o"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/pugixml/vtkpugixml && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/pugixml.dir/src/pugixml.cpp.o -c /home/yyl/SEMO/lib/VTK/ThirdParty/pugixml/vtkpugixml/src/pugixml.cpp

ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/src/pugixml.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pugixml.dir/src/pugixml.cpp.i"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/pugixml/vtkpugixml && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yyl/SEMO/lib/VTK/ThirdParty/pugixml/vtkpugixml/src/pugixml.cpp > CMakeFiles/pugixml.dir/src/pugixml.cpp.i

ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/src/pugixml.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pugixml.dir/src/pugixml.cpp.s"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/pugixml/vtkpugixml && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yyl/SEMO/lib/VTK/ThirdParty/pugixml/vtkpugixml/src/pugixml.cpp -o CMakeFiles/pugixml.dir/src/pugixml.cpp.s

# Object files for target pugixml
pugixml_OBJECTS = \
"CMakeFiles/pugixml.dir/src/pugixml.cpp.o"

# External object files for target pugixml
pugixml_EXTERNAL_OBJECTS =

lib/libvtkpugixml-9.0.so.9.0.1: ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/src/pugixml.cpp.o
lib/libvtkpugixml-9.0.so.9.0.1: ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/build.make
lib/libvtkpugixml-9.0.so.9.0.1: ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/yyl/SEMO/lib/VTK/VTK-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library ../../../lib/libvtkpugixml-9.0.so"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/pugixml/vtkpugixml && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pugixml.dir/link.txt --verbose=$(VERBOSE)
	cd /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/pugixml/vtkpugixml && $(CMAKE_COMMAND) -E cmake_symlink_library ../../../lib/libvtkpugixml-9.0.so.9.0.1 ../../../lib/libvtkpugixml-9.0.so.1 ../../../lib/libvtkpugixml-9.0.so

lib/libvtkpugixml-9.0.so.1: lib/libvtkpugixml-9.0.so.9.0.1
	@$(CMAKE_COMMAND) -E touch_nocreate lib/libvtkpugixml-9.0.so.1

lib/libvtkpugixml-9.0.so: lib/libvtkpugixml-9.0.so.9.0.1
	@$(CMAKE_COMMAND) -E touch_nocreate lib/libvtkpugixml-9.0.so

# Rule to build all files generated by this target.
ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/build: lib/libvtkpugixml-9.0.so

.PHONY : ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/build

ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/clean:
	cd /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/pugixml/vtkpugixml && $(CMAKE_COMMAND) -P CMakeFiles/pugixml.dir/cmake_clean.cmake
.PHONY : ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/clean

ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/depend:
	cd /home/yyl/SEMO/lib/VTK/VTK-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yyl/SEMO/lib/VTK /home/yyl/SEMO/lib/VTK/ThirdParty/pugixml/vtkpugixml /home/yyl/SEMO/lib/VTK/VTK-build /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/pugixml/vtkpugixml /home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ThirdParty/pugixml/vtkpugixml/CMakeFiles/pugixml.dir/depend

