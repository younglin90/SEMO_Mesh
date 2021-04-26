/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkIOExportGL2PSObjectFactory.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkIOExportGL2PSObjectFactory.h"
#include "vtkVersion.h"

// Include all of the classes we want to create overrides for.
#include "vtkOpenGLGL2PSExporter.h"


vtkStandardNewMacro(vtkIOExportGL2PSObjectFactory);

// Now create the functions to create overrides with.
VTK_CREATE_CREATE_FUNCTION(vtkOpenGLGL2PSExporter)


vtkIOExportGL2PSObjectFactory::vtkIOExportGL2PSObjectFactory()
{
this->RegisterOverride("vtkGL2PSExporter", "vtkOpenGLGL2PSExporter", "Override for VTK::IOExportGL2PS module", 1, vtkObjectFactoryCreatevtkOpenGLGL2PSExporter);

}

const char * vtkIOExportGL2PSObjectFactory::GetVTKSourceVersion()
{
  return VTK_SOURCE_VERSION;
}

void vtkIOExportGL2PSObjectFactory::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

// Registration of object factories.
static unsigned int vtkIOExportGL2PSCount = 0;

VTKIOEXPORTGL2PS_EXPORT void vtkIOExportGL2PS_AutoInit_Construct()
{
  if(++vtkIOExportGL2PSCount == 1)
  {


    vtkIOExportGL2PSObjectFactory* factory = vtkIOExportGL2PSObjectFactory::New();
    if (factory)
    {
      // vtkObjectFactory keeps a reference to the "factory",
      vtkObjectFactory::RegisterFactory(factory);
      factory->Delete();
    }
  }
}
