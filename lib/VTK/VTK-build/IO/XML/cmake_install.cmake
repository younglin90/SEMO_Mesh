# Install script for directory: /home/yyl/SEMO/lib/VTK/IO/XML

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
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkRTXMLPolyDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLCompositeDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLCompositeDataWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLDataObjectWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLDataSetWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLFileReadTester.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLGenericDataObjectReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLHierarchicalBoxDataFileConverter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLHierarchicalBoxDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLHierarchicalBoxDataWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLHierarchicalDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLHyperTreeGridReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLHyperTreeGridWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLImageDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLImageDataWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLMultiBlockDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLMultiBlockDataWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLMultiGroupDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPDataObjectReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPHyperTreeGridReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPImageDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPPolyDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPRectilinearGridReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPStructuredDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPStructuredGridReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPTableReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPUnstructuredDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPUnstructuredGridReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPartitionedDataSetCollectionReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPartitionedDataSetCollectionWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPartitionedDataSetReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPartitionedDataSetWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPolyDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLPolyDataWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLRectilinearGridReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLRectilinearGridWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLStructuredDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLStructuredDataWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLStructuredGridReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLStructuredGridWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLTableReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLTableWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLUniformGridAMRReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLUniformGridAMRWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLUnstructuredDataReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLUnstructuredDataWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLUnstructuredGridReader.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLUnstructuredGridWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLWriter.h"
    "/home/yyl/SEMO/lib/VTK/IO/XML/vtkXMLWriterC.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/IO/XML/vtkIOXMLModule.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/vtk/hierarchy/VTK" TYPE FILE RENAME "vtkIOXML-hierarchy.txt" FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkIOXML-hierarchy.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkIOXML-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkIOXML-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkIOXML-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkIOXML-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkIOXML-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkIOXML-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkIOXML-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkIOXML-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkIOXML-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkIOXML-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkIOXML-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkIOXML-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkIOXML-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkIOXML-9.0.so")
    endif()
  endif()
endif()

