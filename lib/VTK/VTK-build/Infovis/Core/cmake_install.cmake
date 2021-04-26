# Install script for directory: /home/yyl/SEMO/lib/VTK/Infovis/Core

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
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkAddMembershipArray.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkAdjacencyMatrixToEdgeTable.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkArrayNorm.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkArrayToTable.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkCollapseGraph.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkCollapseVerticesByArray.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkContinuousScatterplot.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkDataObjectToTable.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkDotProductSimilarity.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkEdgeCenters.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkExpandSelectedGraph.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkExtractSelectedGraph.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkExtractSelectedTree.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkGenerateIndexArray.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkGraphHierarchicalBundleEdges.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkGroupLeafVertices.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkKCoreDecomposition.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkMergeColumns.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkMergeGraphs.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkMergeTables.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkMutableGraphHelper.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkNetworkHierarchy.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkPipelineGraphSource.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkPruneTreeFilter.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkRandomGraphSource.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkReduceTable.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkRemoveHiddenData.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkRemoveIsolatedVertices.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkSparseArrayToTable.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkStreamGraph.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkStringToCategory.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkStringToNumeric.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkTableToArray.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkTableToGraph.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkTableToSparseArray.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkTableToTreeFilter.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkThresholdGraph.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkThresholdTable.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkTransferAttributes.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkTransposeMatrix.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkTreeDifferenceFilter.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkTreeFieldAggregator.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkTreeLevelsFilter.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkVertexDegree.h"
    "/home/yyl/SEMO/lib/VTK/Infovis/Core/vtkWordCloud.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Infovis/Core/vtkInfovisCoreModule.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/vtk/hierarchy/VTK" TYPE FILE RENAME "vtkInfovisCore-hierarchy.txt" FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkInfovisCore-hierarchy.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInfovisCore-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInfovisCore-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkInfovisCore-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkInfovisCore-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInfovisCore-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInfovisCore-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInfovisCore-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInfovisCore-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInfovisCore-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkInfovisCore-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInfovisCore-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInfovisCore-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInfovisCore-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInfovisCore-9.0.so")
    endif()
  endif()
endif()

