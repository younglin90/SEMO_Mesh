
#ifndef VTKFILTERSFLOWPATHS_EXPORT_H
#define VTKFILTERSFLOWPATHS_EXPORT_H

#ifdef VTKFILTERSFLOWPATHS_STATIC_DEFINE
#  define VTKFILTERSFLOWPATHS_EXPORT
#  define VTKFILTERSFLOWPATHS_NO_EXPORT
#else
#  ifndef VTKFILTERSFLOWPATHS_EXPORT
#    ifdef FiltersFlowPaths_EXPORTS
        /* We are building this library */
#      define VTKFILTERSFLOWPATHS_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define VTKFILTERSFLOWPATHS_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef VTKFILTERSFLOWPATHS_NO_EXPORT
#    define VTKFILTERSFLOWPATHS_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef VTKFILTERSFLOWPATHS_DEPRECATED
#  define VTKFILTERSFLOWPATHS_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef VTKFILTERSFLOWPATHS_DEPRECATED_EXPORT
#  define VTKFILTERSFLOWPATHS_DEPRECATED_EXPORT VTKFILTERSFLOWPATHS_EXPORT VTKFILTERSFLOWPATHS_DEPRECATED
#endif

#ifndef VTKFILTERSFLOWPATHS_DEPRECATED_NO_EXPORT
#  define VTKFILTERSFLOWPATHS_DEPRECATED_NO_EXPORT VTKFILTERSFLOWPATHS_NO_EXPORT VTKFILTERSFLOWPATHS_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef VTKFILTERSFLOWPATHS_NO_DEPRECATED
#    define VTKFILTERSFLOWPATHS_NO_DEPRECATED
#  endif
#endif

/* AutoInit implementations. */
#ifdef vtkFiltersFlowPaths_AUTOINIT_INCLUDE
#include vtkFiltersFlowPaths_AUTOINIT_INCLUDE
#endif
#ifdef vtkFiltersFlowPaths_AUTOINIT
#include "vtkAutoInit.h"
VTK_MODULE_AUTOINIT(vtkFiltersFlowPaths)
#endif

#endif /* VTKFILTERSFLOWPATHS_EXPORT_H */
