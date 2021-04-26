/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRenderingUIObjectFactory.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkRenderingUIObjectFactory_h
#define vtkRenderingUIObjectFactory_h

#include "vtkRenderingUIModule.h" // For export macro
#include "vtkObjectFactory.h"

class VTKRENDERINGUI_EXPORT vtkRenderingUIObjectFactory : public vtkObjectFactory
{
public:
  static vtkRenderingUIObjectFactory * New();
  vtkTypeMacro(vtkRenderingUIObjectFactory, vtkObjectFactory);

  const char * GetDescription() override { return "vtkRenderingUI factory overrides."; }

  const char * GetVTKSourceVersion() override;

  void PrintSelf(ostream &os, vtkIndent indent) override;

protected:
  vtkRenderingUIObjectFactory();

private:
  vtkRenderingUIObjectFactory(const vtkRenderingUIObjectFactory&) = delete;
  void operator=(const vtkRenderingUIObjectFactory&) = delete;
};

#endif // vtkRenderingUIObjectFactory_h
