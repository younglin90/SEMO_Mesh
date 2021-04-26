# Install script for directory: /home/yyl/SEMO/lib/VTK/Filters/Parallel

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
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkBlockDistribution.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkAdaptiveTemporalInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkAggregateDataSetFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkAngularPeriodicFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkCollectGraph.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkCollectPolyData.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkCollectTable.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkCutMaterial.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkDistributedDataFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkDuplicatePolyData.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkExtractCTHPart.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkExtractPolyDataPiece.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkExtractUnstructuredGridPiece.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkExtractUserDefinedPiece.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkHyperTreeGridGhostCellsGenerator.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkIntegrateAttributes.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPCellDataToPointData.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPExtractDataArraysOverTime.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPExtractExodusGlobalTemporalVariables.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPExtractSelectedArraysOverTime.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPeriodicFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPKdTree.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPLinearExtrusionFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPMaskPoints.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPMergeArrays.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPOutlineCornerFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPOutlineFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPOutlineFilterInternals.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPPolyDataNormals.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPProbeFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPProjectSphereFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPReflectionFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPResampleFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPSphereSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPYoungsMaterialInterface.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPassThroughFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPieceRequestFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPieceScalars.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPipelineSize.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkProcessIdScalars.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkPTextureMapToSphere.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkRectilinearGridOutlineFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkRemoveGhosts.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkTransmitPolyDataPiece.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkTransmitStructuredDataPiece.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkTransmitRectilinearGridPiece.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkTransmitStructuredGridPiece.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkTransmitUnstructuredGridPiece.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Parallel/vtkUnstructuredGridGhostCellsGenerator.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Filters/Parallel/vtkFiltersParallelModule.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/vtk/hierarchy/VTK" TYPE FILE RENAME "vtkFiltersParallel-hierarchy.txt" FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkFiltersParallel-hierarchy.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersParallel-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersParallel-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersParallel-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersParallel-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersParallel-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersParallel-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersParallel-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersParallel-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersParallel-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersParallel-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersParallel-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersParallel-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersParallel-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersParallel-9.0.so")
    endif()
  endif()
endif()

