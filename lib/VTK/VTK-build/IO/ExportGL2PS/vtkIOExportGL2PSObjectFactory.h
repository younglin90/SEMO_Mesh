/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkIOExportGL2PSObjectFactory.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkIOExportGL2PSObjectFactory_h
#define vtkIOExportGL2PSObjectFactory_h

#include "vtkIOExportGL2PSModule.h" // For export macro
#include "vtkObjectFactory.h"

class VTKIOEXPORTGL2PS_EXPORT vtkIOExportGL2PSObjectFactory : public vtkObjectFactory
{
public:
  static vtkIOExportGL2PSObjectFactory * New();
  vtkTypeMacro(vtkIOExportGL2PSObjectFactory, vtkObjectFactory);

  const char * GetDescription() override { return "vtkIOExportGL2PS factory overrides."; }

  const char * GetVTKSourceVersion() override;

  void PrintSelf(ostream &os, vtkIndent indent) override;

protected:
  vtkIOExportGL2PSObjectFactory();

private:
  vtkIOExportGL2PSObjectFactory(const vtkIOExportGL2PSObjectFactory&) = delete;
  void operator=(const vtkIOExportGL2PSObjectFactory&) = delete;
};

#endif // vtkIOExportGL2PSObjectFactory_h
