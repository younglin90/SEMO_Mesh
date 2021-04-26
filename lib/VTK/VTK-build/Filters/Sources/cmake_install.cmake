# Install script for directory: /home/yyl/SEMO/lib/VTK/Filters/Sources

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
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkArcSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkArrowSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkButtonSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkCapsuleSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkCellTypeSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkConeSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkCubeSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkCylinderSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkDiagonalMatrixSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkDiskSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkEllipseArcSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkEllipticalButtonSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkFrustumSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkGlyphSource2D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkGraphToPolyData.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkHyperTreeGridSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkLineSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkOutlineCornerFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkOutlineCornerSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkOutlineSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkParametricFunctionSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkPlaneSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkPlatonicSolidSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkPointSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkPolyLineSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkPolyPointSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkProgrammableDataObjectSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkProgrammableSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkRandomHyperTreeGridSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkRectangularButtonSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkRegularPolygonSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkSelectionSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkSphereSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkSuperquadricSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkTessellatedBoxSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkTextSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkTexturedSphereSource.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Sources/vtkUniformHyperTreeGridSource.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Filters/Sources/vtkFiltersSourcesModule.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/vtk/hierarchy/VTK" TYPE FILE RENAME "vtkFiltersSources-hierarchy.txt" FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkFiltersSources-hierarchy.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersSources-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersSources-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersSources-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersSources-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersSources-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersSources-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersSources-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersSources-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersSources-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersSources-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersSources-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersSources-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersSources-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersSources-9.0.so")
    endif()
  endif()
endif()

