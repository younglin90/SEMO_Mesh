# Install script for directory: /home/yyl/SEMO/lib/VTK/Common/Core

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
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkABI.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayIteratorIncludes.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkAssume.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkAtomicTypeConcepts.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkAutoInit.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkBuffer.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkCollectionRange.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDataArrayAccessor.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDataArrayIteratorMacro.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDataArrayMeta.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDataArrayRange.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDataArrayTemplate.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDataArrayTupleRange_AOS.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDataArrayTupleRange_Generic.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDataArrayValueRange_AOS.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDataArrayValueRange_Generic.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkEventData.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkGenericDataArrayLookupHelper.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkIOStream.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkIOStreamFwd.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationInternals.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkMathUtilities.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkMeta.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkNew.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkRange.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkRangeIterableTraits.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkSetGet.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkSmartPointer.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkSystemIncludes.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkTemplateAliasMacro.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkTestDataArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkType.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkTypeTraits.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkTypedDataArrayIterator.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkVariantCast.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkVariantCreate.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkVariantExtract.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkVariantInlineOperators.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkWeakPointer.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkWin32Header.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkWindows.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkWrappingHints.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkAtomic.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkSMPThreadLocal.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkSMPToolsInternal.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkSMPTools.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkSMPThreadLocalObject.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkArrayDispatchArrayList.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkConfigure.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkMathConfigure.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkToolkits.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkTypeListMacros.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkVersionMacros.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkVTK_USE_SCALED_SOA_ARRAYS.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkTypeInt8Array.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkTypeInt16Array.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkTypeInt32Array.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkTypeInt64Array.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkTypeUInt8Array.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkTypeUInt16Array.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkTypeUInt32Array.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkTypeUInt64Array.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkTypeFloat32Array.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkTypeFloat64Array.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkAbstractArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkAnimationCue.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArchiver.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayCoordinates.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayExtents.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayExtentsList.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayIterator.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayRange.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArraySort.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayWeights.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkBitArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkBitArrayIterator.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkBoxMuellerRandomSequence.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkBreakPoint.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkByteSwap.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkCallbackCommand.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkCharArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkCollection.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkCollectionIterator.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkCommand.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkCommonInformationKeyManager.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkConditionVariable.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkCriticalSection.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDataArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDataArrayCollection.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDataArrayCollectionIterator.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDataArraySelection.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDebugLeaks.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDebugLeaksManager.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDoubleArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDynamicLoader.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkEventForwarderCommand.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkFileOutputWindow.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkFloatArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkFloatingPointExceptions.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkGarbageCollector.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkGarbageCollectorManager.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkGaussianRandomSequence.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkIdList.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkIdListCollection.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkIdTypeArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkIndent.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformation.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationDataObjectKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationDoubleKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationDoubleVectorKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationIdTypeKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationInformationKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationInformationVectorKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationIntegerKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationIntegerPointerKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationIntegerVectorKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationIterator.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationKeyLookup.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationKeyVectorKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationObjectBaseKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationObjectBaseVectorKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationRequestKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationStringKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationStringVectorKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationUnsignedLongKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationVariantKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationVariantVectorKey.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkInformationVector.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkIntArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkLargeInteger.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkLogger.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkLongArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkLongLongArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkLookupTable.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkMath.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkMersenneTwister.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkMinimalStandardRandomSequence.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkMultiThreader.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkMutexLock.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkOStrStreamWrapper.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkOStreamWrapper.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkObject.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkObjectBase.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkObjectFactory.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkObjectFactoryCollection.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkOldStyleCallbackCommand.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkOutputWindow.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkOverrideInformation.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkOverrideInformationCollection.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkPoints.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkPoints2D.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkPriorityQueue.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkRandomPool.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkRandomSequence.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkReferenceCount.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkScalarsToColors.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkShortArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkSignedCharArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkSimpleCriticalSection.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkSmartPointerBase.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkSortDataArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkStdString.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkStringArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkStringOutputWindow.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkTimePointUtility.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkTimeStamp.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkUnicodeString.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkUnicodeStringArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkUnsignedCharArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkUnsignedIntArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkUnsignedLongArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkUnsignedLongLongArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkUnsignedShortArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkVariant.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkVariantArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkVersion.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkVoidArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkWeakPointerBase.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkWeakReference.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkWindow.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkXMLFileOutputWindow.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkAOSDataArrayTemplate.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayDispatch.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayInterpolate.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayIteratorTemplate.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayPrint.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDenseArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkGenericDataArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkMappedDataArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkSOADataArrayTemplate.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkSparseArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkTypedArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkTypedDataArray.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkTypeList.h"
    "/home/yyl/SEMO/lib/VTK/VTK-build/Common/Core/vtkCommonCoreModule.h"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayIteratorTemplateImplicit.txx"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkAOSDataArrayTemplate.txx"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayDispatch.txx"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayInterpolate.txx"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayIteratorTemplate.txx"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkArrayPrint.txx"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkDenseArray.txx"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkGenericDataArray.txx"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkMappedDataArray.txx"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkSOADataArrayTemplate.txx"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkSparseArray.txx"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkTypedArray.txx"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkTypedDataArray.txx"
    "/home/yyl/SEMO/lib/VTK/Common/Core/vtkTypeList.txx"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/vtk/hierarchy/VTK" TYPE FILE RENAME "vtkCommonCore-hierarchy.txt" FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/vtk/hierarchy/VTK/vtkCommonCore-hierarchy.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkCommonCore-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkCommonCore-9.0.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkCommonCore-9.0.so.9.0.1"
    "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkCommonCore-9.0.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkCommonCore-9.0.so.9.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkCommonCore-9.0.so.1"
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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkCommonCore-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkCommonCore-9.0.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkCommonCore-9.0.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/yyl/SEMO/lib/VTK/VTK-build/lib/libvtkCommonCore-9.0.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkCommonCore-9.0.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkCommonCore-9.0.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkCommonCore-9.0.so"
         OLD_RPATH "/home/yyl/SEMO/lib/VTK/VTK-build/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libvtkCommonCore-9.0.so")
    endif()
  endif()
endif()

