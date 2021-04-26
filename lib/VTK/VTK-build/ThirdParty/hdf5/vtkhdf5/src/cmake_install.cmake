# Install script for directory: /home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/vtk-9.0/vtkhdf5/src" TYPE FILE FILES
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/hdf5.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5api_adpt.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5public.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/vtk_hdf5_mangle.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Apublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5ACpublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Cpublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Dpublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Epublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Fpublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5FDcore.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5FDdirect.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5FDfamily.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5FDhdfs.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5FDlog.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5FDmpi.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5FDmpio.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5FDmulti.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5FDpublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5FDros3.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5FDs3comms.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5FDsec2.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5FDstdio.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5FDwindows.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Gpublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Ipublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Lpublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5MMpublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Opublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Ppublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5PLextern.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5PLpublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Rpublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Spublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Tpublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Zpublic.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Edefin.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Einit.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Epubgen.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5Eterm.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5version.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/hdf5/vtkhdf5/src/H5overflow.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkhdf5-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkhdf5-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkhdf5-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkhdf5-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkhdf5-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkhdf5-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkhdf5-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkhdf5-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkhdf5-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkhdf5-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkhdf5-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkhdf5-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkhdf5-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkhdf5-9.0.so")
    endif()
  endif()
endif()

