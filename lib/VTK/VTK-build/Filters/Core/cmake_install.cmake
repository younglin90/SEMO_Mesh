# Install script for directory: /home/yyl/SEMO/lib/VTK/Filters/Core

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
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtk3DLinearGridPlaneCutter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkAppendArcLength.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkAppendCompositeDataLeaves.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkAppendDataSets.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkAppendFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkAppendPolyData.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkAppendSelection.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkArrayCalculator.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkAssignAttribute.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkAttributeDataToFieldDataFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkBinCellDataFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkCellDataToPointData.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkCellCenters.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkCenterOfMass.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkCleanPolyData.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkClipPolyData.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkCompositeCutter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkCompositeDataProbeFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkConnectivityFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkContour3DLinearGrid.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkContourFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkContourGrid.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkContourHelper.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtk3DLinearGridCrinkleExtractor.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkCutter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkDataObjectGenerator.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkDataObjectToDataSetFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkDataSetEdgeSubdivisionCriterion.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkDataSetToDataObjectFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkDecimatePolylineFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkDecimatePro.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkDelaunay2D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkDelaunay3D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkEdgeSubdivisionCriterion.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkElevationFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkExecutionTimer.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkExplicitStructuredGridCrop.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkExplicitStructuredGridToUnstructuredGrid.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkFeatureEdges.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkFieldDataToAttributeDataFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkFlyingEdges2D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkFlyingEdges3D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkFlyingEdgesPlaneCutter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkGlyph2D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkGlyph3D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkGridSynchronizedTemplates3D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkHedgeHog.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkHull.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkIdFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkImageDataToExplicitStructuredGrid.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkImageAppend.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkImplicitPolyDataDistance.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkImplicitProjectOnPlaneDistance.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkMarchingCubes.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkMarchingSquares.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkMaskFields.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkMaskPoints.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkMaskPolyData.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkMassProperties.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkMergeDataObjectFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkMergeFields.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkMergeFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkMoleculeAppend.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkMultiObjectMassProperties.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkPlaneCutter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkPointDataToCellData.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkPolyDataConnectivityFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkPolyDataNormals.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkPolyDataTangents.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkProbeFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkQuadricClustering.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkQuadricDecimation.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkRearrangeFields.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkRectilinearSynchronizedTemplates.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkRemoveDuplicatePolys.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkResampleToImage.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkResampleWithDataSet.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkReverseSense.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkSimpleElevationFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkSmoothPolyDataFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkSphereTreeFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkStaticCleanPolyData.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkStreamerBase.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkStreamingTessellator.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkStripper.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkStructuredGridAppend.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkStructuredGridOutlineFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkSynchronizedTemplates2D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkSynchronizedTemplates3D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkSynchronizedTemplatesCutter3D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkTensorGlyph.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkThreshold.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkThresholdPoints.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkTransposeTable.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkTriangleFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkTriangleMeshPointNormals.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkTubeFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkUnstructuredGridQuadricDecimation.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkUnstructuredGridToExplicitStructuredGrid.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkVectorDot.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkVectorNorm.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkVoronoi2D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/Core/vtkWindowedSincPolyDataFilter.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Filters/Core/vtkFiltersCoreModule.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/vtk/hierarchy/VTK" TYPE FILE RENAME "vtkFiltersCore-hierarchy.txt" FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkFiltersCore-hierarchy.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersCore-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersCore-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersCore-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersCore-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersCore-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersCore-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersCore-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersCore-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersCore-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersCore-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersCore-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersCore-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersCore-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersCore-9.0.so")
    endif()
  endif()
endif()

