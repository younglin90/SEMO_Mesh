# Install script for directory: /home/yyl/SEMO/lib/VTK/Rendering/Core

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
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkGPUInfoListArray.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkNoise200x200.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkPythagoreanQuadruples.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkRayCastStructures.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkRenderingCoreEnums.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Rendering/Core/vtkTDxConfigure.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkTDxMotionEventInfo.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkAbstractMapper.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkAbstractMapper3D.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkAbstractPicker.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkAbstractVolumeMapper.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkActor.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkActor2D.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkActor2DCollection.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkActorCollection.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkAssembly.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkAvatar.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkBackgroundColorMonitor.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkBillboardTextActor3D.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkCIEDE2000.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkCamera.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkCameraActor.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkCameraInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkCellCenterDepthSort.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkColorTransferFunction.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkCompositeDataDisplayAttributes.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkCompositeDataDisplayAttributesLegacy.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkCompositePolyDataMapper.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkCoordinate.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkCuller.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkCullerCollection.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkDataSetMapper.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkDiscretizableColorTransferFunction.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkDistanceToCamera.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkFXAAOptions.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkFlagpoleLabel.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkFollower.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkFrameBufferObjectBase.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkFrustumCoverageCuller.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkGPUInfo.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkGPUInfoList.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkGenericVertexAttributeMapping.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkGlyph3DMapper.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkGraphMapper.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkGraphToGlyphs.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkGraphicsFactory.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkHardwareSelector.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkHardwareWindow.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkHierarchicalPolyDataMapper.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkImageActor.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkImageMapper.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkImageMapper3D.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkImageProperty.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkImageSlice.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkImageSliceMapper.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkInteractorEventRecorder.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkInteractorObserver.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkLabeledContourMapper.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkLight.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkLightActor.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkLightCollection.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkLightKit.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkLogLookupTable.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkLookupTableWithEnabling.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkMapArrayValues.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkMapper.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkMapper2D.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkMapperCollection.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkObserverMediator.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkPointGaussianMapper.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkPolyDataMapper.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkPolyDataMapper2D.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkProp.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkProp3D.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkProp3DCollection.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkProp3DFollower.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkPropAssembly.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkPropCollection.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkProperty.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkProperty2D.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkRenderPass.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkRenderState.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkRenderTimerLog.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkRenderWindow.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkRenderWindowCollection.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkRenderWindowInteractor.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkRenderWindowInteractor3D.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkRenderer.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkRendererCollection.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkRendererDelegate.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkRendererSource.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkSelectVisiblePoints.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkShaderProperty.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkSkybox.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkStereoCompositor.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkTextActor.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkTextActor3D.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkTexture.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkTexturedActor2D.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkTransformCoordinateSystems.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkTransformInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkTupleInterpolator.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkUniforms.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkViewDependentErrorMetric.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkViewport.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkVisibilitySort.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkVolume.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkVolumeCollection.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkVolumeProperty.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkWindowLevelLookupTable.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkWindowToImageFilter.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkAssemblyNode.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkAssemblyPath.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkAssemblyPaths.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkAreaPicker.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkPicker.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkAbstractPropPicker.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkLODProp3D.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkPropPicker.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkPickingManager.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkWorldPointPicker.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkCellPicker.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkPointPicker.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkRenderedAreaPicker.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkScenePicker.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkInteractorStyle.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkInteractorStyle3D.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkInteractorStyleSwitchBase.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkTDxInteractorStyle.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkTDxInteractorStyleCamera.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkTDxInteractorStyleSettings.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkStringToImage.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkTextMapper.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkTextProperty.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkTextPropertyCollection.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkTextRenderer.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkAbstractInteractionDevice.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkAbstractRenderDevice.h"
    "/home/yyl/SEMO/lib/VTK/Rendering/Core/vtkRenderWidget.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Rendering/Core/vtkRenderingCoreModule.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/vtk/hierarchy/VTK" TYPE FILE RENAME "vtkRenderingCore-hierarchy.txt" FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkRenderingCore-hierarchy.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkRenderingCore-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkRenderingCore-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkRenderingCore-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkRenderingCore-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkRenderingCore-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkRenderingCore-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkRenderingCore-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkRenderingCore-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkRenderingCore-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkRenderingCore-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkRenderingCore-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkRenderingCore-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkRenderingCore-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkRenderingCore-9.0.so")
    endif()
  endif()
endif()

