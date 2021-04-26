# Install script for directory: /home/yyl/SEMO/lib/VTK/Filters/Geometry

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
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkAbstractGridConnectivity.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkCompositeDataGeometryFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkDataSetGhostGenerator.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkDataSetRegionSurfaceFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkDataSetSurfaceFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkExplicitStructuredGridSurfaceFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkGeometryFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkHierarchicalDataSetGeometryFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkImageDataGeometryFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkImageDataToUniformGrid.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkLinearToQuadraticCellsFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkProjectSphereFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkRectilinearGridGeometryFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkRectilinearGridPartitioner.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkStructuredAMRGridConnectivity.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkStructuredAMRNeighbor.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkStructuredGridConnectivity.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkStructuredGridGeometryFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkStructuredGridGhostDataGenerator.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkStructuredGridPartitioner.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkStructuredNeighbor.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkStructuredPointsGeometryFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkUniformGridGhostDataGenerator.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Geometry/vtkUnstructuredGridGeometryFilter.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Filters/Geometry/vtkFiltersGeometryModule.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/vtk/hierarchy/VTK" TYPE FILE RENAME "vtkFiltersGeometry-hierarchy.txt" FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkFiltersGeometry-hierarchy.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeometry-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeometry-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersGeometry-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersGeometry-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeometry-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeometry-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeometry-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeometry-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeometry-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersGeometry-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeometry-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeometry-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeometry-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeometry-9.0.so")
    endif()
  endif()
endif()

