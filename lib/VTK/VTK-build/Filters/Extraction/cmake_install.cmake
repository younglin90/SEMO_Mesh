# Install script for directory: /home/yyl/SEMO/lib/VTK/Filters/Extraction

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/yyl/SEMO/lib/VTK/VTK-build")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/vtk-9.0" TYPE FILE FILES
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkBlockSelector.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkConvertSelection.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExpandMarkedElements.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractBlock.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractCells.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractCellsByType.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractDataArraysOverTime.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractDataOverTime.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractDataSets.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractEdges.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractExodusGlobalTemporalVariables.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractGeometry.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractGrid.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractLevel.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractPolyDataGeometry.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractRectilinearGrid.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractSelectedArraysOverTime.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractSelectedBlock.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractSelectedIds.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractSelectedLocations.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractSelectedPolyDataIds.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractSelectedRows.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractSelectedThresholds.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractSelection.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractTensorComponents.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractTimeSteps.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractUnstructuredGrid.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractVectorComponents.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkFrustumSelector.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkHierarchicalDataExtractDataSets.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkHierarchicalDataExtractLevel.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkLocationSelector.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkProbeSelectedLocations.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkSelector.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkValueSelector.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Extraction/vtkExtractTemporalFieldData.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Filters/Extraction/vtkFiltersExtractionModule.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/vtk/hierarchy/VTK" TYPE FILE RENAME "vtkFiltersExtraction-hierarchy.txt" FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkFiltersExtraction-hierarchy.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersExtraction-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersExtraction-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersExtraction-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersExtraction-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersExtraction-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersExtraction-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
           NEW_RPATH "")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersExtraction-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersExtraction-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersExtraction-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersExtraction-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersExtraction-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersExtraction-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersExtraction-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersExtraction-9.0.so")
    endif()
  endif()
endif()

