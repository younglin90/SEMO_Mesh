# Install script for directory: /home/yyl/SEMO/lib/VTK/Charts/Core

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
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkAxis.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkAxisExtended.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkCategoryLegend.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkChartBox.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkChart.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkChartHistogram2D.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkChartLegend.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkChartMatrix.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkChartParallelCoordinates.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkChartPie.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkChartXY.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkChartXYZ.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkColorLegend.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkColorTransferControlPointsItem.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkColorTransferFunctionItem.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkCompositeControlPointsItem.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkCompositeTransferFunctionItem.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkContextArea.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkContextPolygon.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkControlPointsItem.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkInteractiveArea.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkLookupTableItem.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPiecewiseControlPointsItem.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPiecewiseFunctionItem.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPiecewisePointHandleItem.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlot3D.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlotArea.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlotBag.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlotBar.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlotBox.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlot.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlotFunctionalBag.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlotGrid.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlotHistogram2D.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlotLine3D.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlotLine.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlotParallelCoordinates.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlotPie.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlotPoints3D.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlotPoints.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlotStacked.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkPlotSurface.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkRangeHandlesItem.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkScalarsToColorsItem.h"
    "/home/yyl/SEMO/lib/VTK/Charts/Core/vtkScatterPlotMatrix.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Charts/Core/vtkChartsCoreModule.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/vtk/hierarchy/VTK" TYPE FILE RENAME "vtkChartsCore-hierarchy.txt" FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkChartsCore-hierarchy.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkChartsCore-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkChartsCore-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkChartsCore-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkChartsCore-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkChartsCore-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkChartsCore-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkChartsCore-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkChartsCore-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkChartsCore-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkChartsCore-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkChartsCore-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkChartsCore-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkChartsCore-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkChartsCore-9.0.so")
    endif()
  endif()
endif()

