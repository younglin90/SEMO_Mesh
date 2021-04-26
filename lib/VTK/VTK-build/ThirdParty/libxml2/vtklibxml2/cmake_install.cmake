# Install script for directory: /home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/vtk-9.0/vtklibxml2/include/libxml" TYPE FILE FILES
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/c14n.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/catalog.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/chvalid.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/debugXML.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/dict.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/DOCBparser.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/encoding.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/entities.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/globals.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/hash.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/HTMLparser.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/HTMLtree.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/list.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/nanoftp.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/nanohttp.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/parser.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/parserInternals.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/pattern.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/relaxng.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/SAX.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/SAX2.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/schemasInternals.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/schematron.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/threads.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/tree.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/uri.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/valid.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/vtk_libxml2_mangle.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xinclude.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xlink.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlautomata.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlerror.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlexports.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlIO.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlmemory.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlmodule.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlreader.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlregexp.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlsave.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlschemas.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlschemastypes.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlstring.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlunicode.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlwin32version.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlwriter.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xpath.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xpathInternals.h"
    "/home/yyl/SEMO/lib/VTK/ThirdParty/libxml2/vtklibxml2/include/libxml/xpointer.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/ThirdParty/libxml2/vtklibxml2/include/libxml/xmlversion.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibxml2-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibxml2-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtklibxml2-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtklibxml2-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibxml2-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibxml2-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibxml2-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibxml2-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibxml2-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtklibxml2-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibxml2-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibxml2-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibxml2-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtklibxml2-9.0.so")
    endif()
  endif()
endif()

