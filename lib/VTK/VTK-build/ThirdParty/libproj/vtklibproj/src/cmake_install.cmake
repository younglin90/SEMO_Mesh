# Install script for directory: /home/yyl/SEMO/lib/VTK/ThirdParty/libproj/vtklibproj/src

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/vtk-9.0/vtklibproj/src" TYPE FILE FILES
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libproj/vtklibproj/src/projects.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libproj/vtklibproj/src/proj_api.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libproj/vtklibproj/src/geodesic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libproj/vtklibproj/src/vtk_libproj_mangle.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/libproj/vtklibproj/src/vtklibproj_export.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibproj-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibproj-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtklibproj-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtklibproj-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibproj-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibproj-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibproj-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibproj-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibproj-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtklibproj-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibproj-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibproj-9.0.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibproj-9.0.so")
    endif()
  endif()
endif()
