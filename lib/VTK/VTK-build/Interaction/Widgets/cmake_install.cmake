# Install script for directory: /home/yyl/SEMO/lib/VTK/Interaction/Widgets

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
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtk3DWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkAbstractPolygonalHandleRepresentation3D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkAbstractWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkAffineRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkAffineRepresentation2D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkAffineWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkAngleRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkAngleRepresentation2D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkAngleRepresentation3D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkAngleWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkAxesTransformRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkAxesTransformWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkBalloonRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkBalloonWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkBezierContourLineInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkBiDimensionalRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkBiDimensionalRepresentation2D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkBiDimensionalWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkBorderRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkBorderWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkBoundedPlanePointPlacer.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkBoxRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkBoxWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkBoxWidget2.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkBrokenLineWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkButtonRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkButtonWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkCameraRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkCameraWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkCaptionRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkCaptionWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkCellCentersPointPlacer.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkCenteredSliderRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkCenteredSliderWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkCheckerboardRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkCheckerboardWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkClosedSurfacePointPlacer.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkConstrainedPointHandleRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkContinuousValueWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkContinuousValueWidgetRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkContourLineInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkContourRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkContourWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkCurveRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkDijkstraImageContourLineInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkDistanceRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkDistanceRepresentation2D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkDistanceRepresentation3D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkDistanceWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkEllipsoidTensorProbeRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkEvent.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkFinitePlaneRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkFinitePlaneWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkFixedSizeHandleRepresentation3D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkFocalPlaneContourRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkFocalPlanePointPlacer.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkHandleRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkHandleWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkHoverWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkImageActorPointPlacer.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkImageCroppingRegionsWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkImageOrthoPlanes.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkImagePlaneWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkImageTracerWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkImplicitCylinderRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkImplicitCylinderWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkImplicitPlaneRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkImplicitPlaneWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkImplicitPlaneWidget2.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkLightRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkLightWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkLinearContourLineInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkLineRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkLineWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkLineWidget2.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkLogoRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkLogoWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkMeasurementCubeHandleRepresentation3D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkOrientationMarkerWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkOrientedGlyphContourRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkOrientedGlyphFocalPlaneContourRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkOrientedPolygonalHandleRepresentation3D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkParallelopipedRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkParallelopipedWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkPlaneWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkPlaybackRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkPlaybackWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkPointHandleRepresentation2D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkPointHandleRepresentation3D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkPointPlacer.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkPointWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkPolyDataContourLineInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkPolyDataPointPlacer.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkPolyDataSourceWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkPolygonalHandleRepresentation3D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkPolygonalSurfaceContourLineInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkPolygonalSurfacePointPlacer.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkPolyLineRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkPolyLineWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkProgressBarRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkProgressBarWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkProp3DButtonRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkRectilinearWipeRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkRectilinearWipeWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkResliceCursor.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkResliceCursorActor.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkResliceCursorLineRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkResliceCursorPicker.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkResliceCursorPolyDataAlgorithm.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkResliceCursorRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkResliceCursorThickLineRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkResliceCursorWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkScalarBarRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkScalarBarWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkSeedRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkSeedWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkSliderRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkSliderRepresentation2D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkSliderRepresentation3D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkSliderWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkSphereHandleRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkSphereRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkSphereWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkSphereWidget2.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkSplineRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkSplineWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkSplineWidget2.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkTensorProbeRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkTensorProbeWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkTerrainContourLineInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkTerrainDataPointPlacer.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkTextRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkTexturedButtonRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkTexturedButtonRepresentation2D.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkTextWidget.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkWidgetCallbackMapper.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkWidgetEvent.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkWidgetEventTranslator.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkWidgetRepresentation.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkWidgetSet.h"
    "/home/yyl/SEMO/lib/VTK/Interaction/Widgets/vtkXYPlotWidget.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Interaction/Widgets/vtkInteractionWidgetsModule.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/vtk/hierarchy/VTK" TYPE FILE RENAME "vtkInteractionWidgets-hierarchy.txt" FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkInteractionWidgets-hierarchy.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInteractionWidgets-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInteractionWidgets-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkInteractionWidgets-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkInteractionWidgets-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInteractionWidgets-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInteractionWidgets-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInteractionWidgets-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInteractionWidgets-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInteractionWidgets-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkInteractionWidgets-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInteractionWidgets-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInteractionWidgets-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInteractionWidgets-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkInteractionWidgets-9.0.so")
    endif()
  endif()
endif()

