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

# Utility rule file for vtkIOImport-hierarchy.

# Include the progress variables for this target.
include IO/Import/CMakeFiles/vtkIOImport-hierarchy.dir/progress.make

IO/Import/CMakeFiles/vtkIOImport-hierarchy: lib/vtk/hierarchy/VTK/vtkIOImport-hierarchy.txt
IO/Import/CMakeFiles/vtkIOImport-hierarchy: bin/vtkWrapHierarchy-9.0


lib/vtk/hierarchy/VTK/vtkIOImport-hierarchy.txt: ../IO/Import/vtk3DS.h
lib/vtk/hierarchy/VTK/vtkIOImport-hierarchy.txt: ../IO/Import/vtk3DSImporter.h
lib/vtk/hierarchy/VTK/vtkIOImport-hierarchy.txt: ../IO/Import/vtkGLTFImporter.h
lib/vtk/hierarchy/VTK/vtkIOImport-hierarchy.txt: ../IO/Import/vtkImporter.h
lib/vtk/hierarchy/VTK/vtkIOImport-hierarchy.txt: ../IO/Import/vtkVRMLImporter.h
lib/vtk/hierarchy/VTK/vtkIOImport-hierarchy.txt: ../IO/Import/vtkOBJImporter.h
lib/vtk/hierarchy/VTK/vtkIOImport-hierarchy.txt: ../IO/Import/vtkOBJImporterInternals.h
lib/vtk/hierarchy/VTK/vtkIOImport-hierarchy.txt: IO/Import/vtkIOImportModule.h
lib/vtk/hierarchy/VTK/vtkIOImport-hierarchy.txt: IO/Import/CMakeFiles/vtkIOImport-hierarchy.Debug.args
lib/vtk/hierarchy/VTK/vtkIOImport-hierarchy.txt: IO/Import/CMakeFiles/vtkIOImport-hierarchy.data
lib/vtk/hierarchy/VTK/vtkIOImport-hierarchy.txt: IO/Import/CMakeFiles/vtkIOImport-hierarchy.depends.args
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/yyl/SEMO/lib/VTK/VTK-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating the wrap hierarchy for VTK::IOImport"
	cd /home/yyl/SEMO/lib/VTK/VTK-build/IO/Import && ../../bin/vtkWrapHierarchy-9.0 @/home/yyl/SEMO/lib/VTK/VTK-build/IO/Import/CMakeFiles/vtkIOImport-hierarchy.Debug.args -o /home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkIOImport-hierarchy.txt /home/yyl/SEMO/lib/VTK/VTK-build/IO/Import/CMakeFiles/vtkIOImport-hierarchy.data @/home/yyl/SEMO/lib/VTK/VTK-build/IO/Import/CMakeFiles/vtkIOImport-hierarchy.depends.args

vtkIOImport-hierarchy: IO/Import/CMakeFiles/vtkIOImport-hierarchy
vtkIOImport-hierarchy: lib/vtk/hierarchy/VTK/vtkIOImport-hierarchy.txt
vtkIOImport-hierarchy: IO/Import/CMakeFiles/vtkIOImport-hierarchy.dir/build.make

.PHONY : vtkIOImport-hierarchy

# Rule to build all files generated by this target.
IO/Import/CMakeFiles/vtkIOImport-hierarchy.dir/build: vtkIOImport-hierarchy

.PHONY : IO/Import/CMakeFiles/vtkIOImport-hierarchy.dir/build

IO/Import/CMakeFiles/vtkIOImport-hierarchy.dir/clean:
	cd /home/yyl/SEMO/lib/VTK/VTK-build/IO/Import && $(CMAKE_COMMAND) -P CMakeFiles/vtkIOImport-hierarchy.dir/cmake_clean.cmake
.PHONY : IO/Import/CMakeFiles/vtkIOImport-hierarchy.dir/clean

IO/Import/CMakeFiles/vtkIOImport-hierarchy.dir/depend:
	cd /home/yyl/SEMO/lib/VTK/VTK-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yyl/SEMO/lib/VTK /home/yyl/SEMO/lib/VTK/IO/Import /home/yyl/SEMO/lib/VTK/VTK-build /home/yyl/SEMO/lib/VTK/VTK-build/IO/Import /home/yyl/SEMO/lib/VTK/VTK-build/IO/Import/CMakeFiles/vtkIOImport-hierarchy.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : IO/Import/CMakeFiles/vtkIOImport-hierarchy.dir/depend

