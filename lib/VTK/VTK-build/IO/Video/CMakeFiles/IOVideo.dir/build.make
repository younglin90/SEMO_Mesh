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
include IO/Video/CMakeFiles/IOVideo.dir/depend.make

# Include the progress variables for this target.
include IO/Video/CMakeFiles/IOVideo.dir/progress.make

# Include the compile flags for this target's objects.
include IO/Video/CMakeFiles/IOVideo.dir/flags.make

IO/Video/CMakeFiles/IOVideo.dir/vtkVideoSource.cxx.o: IO/Video/CMakeFiles/IOVideo.dir/flags.make
IO/Video/CMakeFiles/IOVideo.dir/vtkVideoSource.cxx.o: ../IO/Video/vtkVideoSource.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yyl/SEMO/lib/VTK/VTK-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object IO/Video/CMakeFiles/IOVideo.dir/vtkVideoSource.cxx.o"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/IO/Video && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/IOVideo.dir/vtkVideoSource.cxx.o -c /home/yyl/SEMO/lib/VTK/IO/Video/vtkVideoSource.cxx

IO/Video/CMakeFiles/IOVideo.dir/vtkVideoSource.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/IOVideo.dir/vtkVideoSource.cxx.i"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/IO/Video && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yyl/SEMO/lib/VTK/IO/Video/vtkVideoSource.cxx > CMakeFiles/IOVideo.dir/vtkVideoSource.cxx.i

IO/Video/CMakeFiles/IOVideo.dir/vtkVideoSource.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/IOVideo.dir/vtkVideoSource.cxx.s"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/IO/Video && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yyl/SEMO/lib/VTK/IO/Video/vtkVideoSource.cxx -o CMakeFiles/IOVideo.dir/vtkVideoSource.cxx.s

# Object files for target IOVideo
IOVideo_OBJECTS = \
"CMakeFiles/IOVideo.dir/vtkVideoSource.cxx.o"

# External object files for target IOVideo
IOVideo_EXTERNAL_OBJECTS =

lib/libvtkIOVideo-9.0.so.9.0.1: IO/Video/CMakeFiles/IOVideo.dir/vtkVideoSource.cxx.o
lib/libvtkIOVideo-9.0.so.9.0.1: IO/Video/CMakeFiles/IOVideo.dir/build.make
lib/libvtkIOVideo-9.0.so.9.0.1: lib/libvtkCommonExecutionModel-9.0.so.9.0.1
lib/libvtkIOVideo-9.0.so.9.0.1: lib/libvtkCommonDataModel-9.0.so.9.0.1
lib/libvtkIOVideo-9.0.so.9.0.1: lib/libvtkCommonSystem-9.0.so.9.0.1
lib/libvtkIOVideo-9.0.so.9.0.1: lib/libvtkCommonTransforms-9.0.so.9.0.1
lib/libvtkIOVideo-9.0.so.9.0.1: lib/libvtkCommonMath-9.0.so.9.0.1
lib/libvtkIOVideo-9.0.so.9.0.1: lib/libvtkCommonCore-9.0.so.9.0.1
lib/libvtkIOVideo-9.0.so.9.0.1: lib/libvtksys-9.0.so.9.0.1
lib/libvtkIOVideo-9.0.so.9.0.1: IO/Video/CMakeFiles/IOVideo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/yyl/SEMO/lib/VTK/VTK-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library ../../lib/libvtkIOVideo-9.0.so"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/IO/Video && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/IOVideo.dir/link.txt --verbose=$(VERBOSE)
	cd /home/yyl/SEMO/lib/VTK/VTK-build/IO/Video && $(CMAKE_COMMAND) -E cmake_symlink_library ../../lib/libvtkIOVideo-9.0.so.9.0.1 ../../lib/libvtkIOVideo-9.0.so.1 ../../lib/libvtkIOVideo-9.0.so

lib/libvtkIOVideo-9.0.so.1: lib/libvtkIOVideo-9.0.so.9.0.1
	@$(CMAKE_COMMAND) -E touch_nocreate lib/libvtkIOVideo-9.0.so.1

lib/libvtkIOVideo-9.0.so: lib/libvtkIOVideo-9.0.so.9.0.1
	@$(CMAKE_COMMAND) -E touch_nocreate lib/libvtkIOVideo-9.0.so

# Rule to build all files generated by this target.
IO/Video/CMakeFiles/IOVideo.dir/build: lib/libvtkIOVideo-9.0.so

.PHONY : IO/Video/CMakeFiles/IOVideo.dir/build

IO/Video/CMakeFiles/IOVideo.dir/clean:
	cd /home/yyl/SEMO/lib/VTK/VTK-build/IO/Video && $(CMAKE_COMMAND) -P CMakeFiles/IOVideo.dir/cmake_clean.cmake
.PHONY : IO/Video/CMakeFiles/IOVideo.dir/clean

IO/Video/CMakeFiles/IOVideo.dir/depend:
	cd /home/yyl/SEMO/lib/VTK/VTK-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yyl/SEMO/lib/VTK /home/yyl/SEMO/lib/VTK/IO/Video /home/yyl/SEMO/lib/VTK/VTK-build /home/yyl/SEMO/lib/VTK/VTK-build/IO/Video /home/yyl/SEMO/lib/VTK/VTK-build/IO/Video/CMakeFiles/IOVideo.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : IO/Video/CMakeFiles/IOVideo.dir/depend

