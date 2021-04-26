# Install script for directory: /home/yyl/SEMO/lib/VTK/Views/Infovis

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
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkApplyColors.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkApplyIcons.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkDendrogramItem.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkGraphItem.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkGraphLayoutView.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkHeatmapItem.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkHierarchicalGraphPipeline.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkHierarchicalGraphView.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkIcicleView.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkInteractorStyleAreaSelectHover.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkInteractorStyleTreeMapHover.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkParallelCoordinatesHistogramRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkParallelCoordinatesRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkParallelCoordinatesView.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkRenderedGraphRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkRenderedHierarchyRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkRenderedRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkRenderedSurfaceRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkRenderedTreeAreaRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkRenderView.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkSCurveSpline.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkTanglegramItem.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkTreeAreaView.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkTreeHeatmapItem.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkTreeMapView.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkTreeRingView.h"
    "/home/yyl/SEMO/lib/VTK/Views/Infovis/vtkViewUpdater.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Views/Infovis/vtkViewsInfovisModule.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/vtk/hierarchy/VTK" TYPE FILE RENAME "vtkViewsInfovis-hierarchy.txt" FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkViewsInfovis-hierarchy.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkViewsInfovis-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkViewsInfovis-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkViewsInfovis-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkViewsInfovis-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkViewsInfovis-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkViewsInfovis-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkViewsInfovis-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkViewsInfovis-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkViewsInfovis-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkViewsInfovis-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkViewsInfovis-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkViewsInfovis-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkViewsInfovis-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkViewsInfovis-9.0.so")
    endif()
  endif()
endif()

