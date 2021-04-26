# Install script for directory: /home/yyl/SEMO/lib/VTK/Imaging/General

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
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageAnisotropicDiffusion2D.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageAnisotropicDiffusion3D.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageCheckerboard.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageCityBlockDistance.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageConvolve.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageCorrelation.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageEuclideanDistance.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageEuclideanToPolar.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageGaussianSmooth.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageGradient.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageGradientMagnitude.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageHybridMedian2D.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageLaplacian.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageMedian3D.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageNormalize.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageRange3D.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageSeparableConvolution.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageSobel2D.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageSobel3D.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageSpatialAlgorithm.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageVariance3D.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkSimpleImageFilterExample.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageSlab.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/General/vtkImageSlabReslice.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Imaging/General/vtkImagingGeneralModule.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/vtk/hierarchy/VTK" TYPE FILE RENAME "vtkImagingGeneral-hierarchy.txt" FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkImagingGeneral-hierarchy.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingGeneral-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingGeneral-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkImagingGeneral-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkImagingGeneral-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingGeneral-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingGeneral-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingGeneral-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingGeneral-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingGeneral-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkImagingGeneral-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingGeneral-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingGeneral-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingGeneral-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingGeneral-9.0.so")
    endif()
  endif()
endif()

