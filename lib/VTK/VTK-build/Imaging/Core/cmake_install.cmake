# Install script for directory: /home/yyl/SEMO/lib/VTK/Imaging/Core

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
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkAbstractImageInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkExtractVOI.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageAppendComponents.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageBlend.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageBSplineCoefficients.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageBSplineInternals.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageBSplineInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageCacheFilter.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageCast.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageChangeInformation.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageClip.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageConstantPad.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageDataStreamer.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageDecomposeFilter.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageDifference.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageExtractComponents.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageFlip.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageIterateFilter.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageMagnify.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageMapToColors.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageMask.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageMirrorPad.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImagePadFilter.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImagePermute.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImagePointDataIterator.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImagePointIterator.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageResample.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageResize.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageReslice.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageResliceToColors.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageShiftScale.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageShrink3D.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageSincInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageStencilAlgorithm.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageStencilData.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageStencilIterator.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageStencilSource.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageThreshold.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageTranslateExtent.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkImageWrapPad.h"
    "/home/yyl/SEMO/lib/VTK/Imaging/Core/vtkRTAnalyticSource.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Imaging/Core/vtkImagingCoreModule.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/vtk/hierarchy/VTK" TYPE FILE RENAME "vtkImagingCore-hierarchy.txt" FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkImagingCore-hierarchy.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingCore-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingCore-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkImagingCore-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkImagingCore-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingCore-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingCore-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingCore-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingCore-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingCore-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkImagingCore-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingCore-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingCore-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingCore-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkImagingCore-9.0.so")
    endif()
  endif()
endif()

