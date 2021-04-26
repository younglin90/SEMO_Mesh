# Install script for directory: /home/yyl/SEMO/lib/VTK/Filters/General

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
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkAnnotationLink.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkAppendLocationAttributes.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkAppendPoints.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkApproximatingSubdivisionFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkAreaContourSpectrumFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkAxes.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkBlankStructuredGrid.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkBlankStructuredGridWithImage.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkBlockIdScalars.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkBooleanOperationPolyDataFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkBoxClipDataSet.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkBrownianPoints.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkCellDerivatives.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkCellTreeLocator.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkCellValidator.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkClipClosedSurface.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkClipConvexPolyData.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkClipDataSet.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkClipVolume.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkCoincidentPoints.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkContourTriangulator.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkCountFaces.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkCountVertices.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkCursor2D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkCursor3D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkCurvatures.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkDataSetGradient.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkDataSetGradientPrecompute.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkDataSetTriangleFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkDateToNumeric.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkDeformPointSet.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkDensifyPolyData.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkDicer.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkDiscreteFlyingEdges2D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkDiscreteFlyingEdges3D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkDiscreteFlyingEdgesClipper2D.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkDiscreteMarchingCubes.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkDistancePolyDataFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkEdgePoints.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkExtractArray.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkExtractSelectedFrustum.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkExtractSelectionBase.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkGradientFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkGraphLayoutFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkGraphToPoints.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkGraphWeightEuclideanDistanceFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkGraphWeightFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkHierarchicalDataLevelFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkHyperStreamline.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkIconGlyphFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkImageDataToPointSet.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkImageMarchingCubes.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkInterpolateDataSetAttributes.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkInterpolatingSubdivisionFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkIntersectionPolyDataFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkLevelIdScalars.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkLinkEdgels.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkLoopBooleanPolyDataFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkMarchingContourFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkMatricizeArray.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkMergeArrays.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkMergeCells.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkMultiBlockDataGroupFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkMultiBlockFromTimeSeriesFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkMultiBlockMergeFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkMultiThreshold.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkNormalizeMatrixVectors.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkOBBDicer.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkOBBTree.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkOverlappingAMRLevelIdScalars.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkPassArrays.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkPassSelectedArrays.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkPassThrough.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkPointConnectivityFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkPolyDataStreamer.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkPolyDataToReebGraphFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkProbePolyhedron.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkQuadraturePointInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkQuadraturePointsGenerator.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkQuadratureSchemeDictionaryGenerator.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkQuantizePolyDataPoints.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkRandomAttributeGenerator.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkRectilinearGridClip.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkRectilinearGridToPointSet.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkRectilinearGridToTetrahedra.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkRecursiveDividingCubes.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkReflectionFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkRotationFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkSampleImplicitFunctionFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkShrinkFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkShrinkPolyData.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkSpatialRepresentationFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkSplineFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkSplitByCellScalarFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkSplitColumnComponents.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkSplitField.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkStructuredGridClip.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkSubPixelPositionEdgels.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkSubdivisionFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkSynchronizeTimeFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkTableBasedClipDataSet.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkTableToPolyData.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkTableToStructuredGrid.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkTemporalPathLineFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkTemporalStatistics.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkTessellatorFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkTimeSourceExample.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkTransformFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkTransformPolyDataFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkUncertaintyTubeFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkVertexGlyphFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkVolumeContourSpectrumFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkVoxelContoursToSurfaceFilter.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkWarpLens.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkWarpScalar.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkWarpTo.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkWarpVector.h"
    "/home/yyl/SEMO/lib/VTK/Filters/General/vtkYoungsMaterialInterface.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Filters/General/vtkFiltersGeneralModule.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/vtk/hierarchy/VTK" TYPE FILE RENAME "vtkFiltersGeneral-hierarchy.txt" FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkFiltersGeneral-hierarchy.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeneral-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeneral-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersGeneral-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersGeneral-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeneral-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeneral-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeneral-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeneral-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeneral-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkFiltersGeneral-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeneral-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeneral-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeneral-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkFiltersGeneral-9.0.so")
    endif()
  endif()
endif()

