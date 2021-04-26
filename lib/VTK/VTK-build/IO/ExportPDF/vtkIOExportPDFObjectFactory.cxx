/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkIOExportPDFObjectFactory.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkIOExportPDFObjectFactory.h"
#include "vtkVersion.h"

// Include all of the classes we want to create overrides for.
#include "vtkPDFContextDevice2D.h"


vtkStandardNewMacro(vtkIOExportPDFObjectFactory);

// Now create the functions to create overrides with.
VTK_CREATE_CREATE_FUNCTION(vtkPDFContextDevice2D)


vtkIOExportPDFObjectFactory::vtkIOExportPDFObjectFactory()
{
this->RegisterOverride("vtkContextDevice2D", "vtkPDFContextDevice2D", "Override for VTK::IOExportPDF module", 1, vtkObjectFactoryCreatevtkPDFContextDevice2D);

}

const char * vtkIOExportPDFObjectFactory::GetVTKSourceVersion()
{
  return VTK_SOURCE_VERSION;
}

void vtkIOExportPDFObjectFactory::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

// Registration of object factories.
static unsigned int vtkIOExportPDFCount = 0;

VTKIOEXPORTPDF_EXPORT void vtkIOExportPDF_AutoInit_Construct()
{
  if(++vtkIOExportPDFCount == 1)
  {


    vtkIOExportPDFObjectFactory* factory = vtkIOExportPDFObjectFactory::New();
    if (factory)
    {
      // vtkObjectFactory keeps a reference to the "factory",
      vtkObjectFactory::RegisterFactory(factory);
      factory->Delete();
    }
  }
}
