# Install script for directory: /home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/vtk-9.0/vtkmetaio" TYPE FILE FILES
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/localMetaConfiguration.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaArray.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaArrow.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaBlob.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaCommand.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaContour.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaDTITube.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaEllipse.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaEvent.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaFEMObject.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaForm.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaGaussian.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaGroup.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaImage.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaImageTypes.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaImageUtils.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaITKUtils.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaLandmark.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaLine.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaMesh.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaObject.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaOutput.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaScene.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaSurface.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaTransform.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaTubeGraph.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaTube.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaTypes.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaUtils.h"
    "/home/yyl/SEMO/lib/VTK/Utilities/MetaIO/vtkmetaio/metaVesselTube.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Utilities/MetaIO/vtkmetaio/metaIOConfig.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkmetaio-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkmetaio-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkmetaio-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkmetaio-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkmetaio-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkmetaio-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkmetaio-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkmetaio-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkmetaio-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkmetaio-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkmetaio-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkmetaio-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkmetaio-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkmetaio-9.0.so")
    endif()
  endif()
endif()

