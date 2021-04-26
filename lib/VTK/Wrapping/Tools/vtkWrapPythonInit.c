#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void CreateInitFile(const char* libName, FILE* fout)
{
  const char* prefix = "";
  const char* dllexp = "VTK_ABI_EXPORT ";

  fprintf(fout, "// Generated by vtkWrapPythonInit in VTK/Wrapping\n");
  fprintf(fout, "#include \"vtkPython.h\"\n");
  fprintf(fout, "#include \"vtkPythonCompatibility.h\"\n");
  fprintf(fout, "#include \"vtkSystemIncludes.h\"\n");
  fprintf(fout,
    "// Handle compiler warning messages, etc.\n"
    "#if defined( _MSC_VER ) && !defined(VTK_DISPLAY_WIN32_WARNINGS)\n"
    "#pragma warning ( disable : 4706 )\n"
    "#endif // Windows Warnings\n\n");

  fprintf(fout, "extern \"C\" { PyObject *real_init%s(const char * /*unused*/); }\n\n", libName);

  fprintf(fout, "#ifdef VTK_PY3K\n");

  fprintf(fout, "extern  \"C\" { %sPyObject *PyInit_%s%s(); }\n\n", dllexp, prefix, libName);

  fprintf(fout, "PyObject *PyInit_%s()\n", libName);
  fprintf(fout, "{\n");
  fprintf(fout, "  return real_init%s(nullptr);\n", libName);
  fprintf(fout, "}\n");

  fprintf(fout, "#else\n");

  fprintf(fout, "extern  \"C\" { %svoid init%s%s(); }\n\n", dllexp, prefix, libName);

  fprintf(fout, "void init%s()\n", libName);
  fprintf(fout, "{\n");
  fprintf(fout, "  real_init%s(nullptr);\n", libName);
  fprintf(fout, "}\n");

  fprintf(fout, "#endif\n");
}

/* warning this code is also in getclasses.cxx under pcmaker */
/* this routine creates the init file */
static void CreateImplFile(const char* libName, const char* importName, int numDepends,
  char** depends, int numFiles, char** files, FILE* fout)
{
  int i;

  const char* dllexp = "VTK_ABI_EXPORT ";

  fprintf(fout, "// Generated by vtkWrapPythonInit in VTK/Wrapping\n");
  fprintf(fout, "#include \"vtkPythonUtil.h\"\n");
  fprintf(fout, "#include \"vtkSystemIncludes.h\"\n");
  fprintf(fout, "#include <cstring>\n");
  fprintf(fout,
    "// Handle compiler warning messages, etc.\n"
    "#if defined( _MSC_VER ) && !defined(VTK_DISPLAY_WIN32_WARNINGS)\n"
    "#pragma warning ( disable : 4706 )\n"
    "#endif // Windows Warnings\n\n");

  for (i = 0; i < numFiles; i++)
  {
    fprintf(fout, "extern \"C\" { PyObject *PyVTKAddFile_%s(PyObject *dict); }\n", files[i]);
  }

  fprintf(fout, "\nstatic PyMethodDef Py%s_Methods[] = {\n", libName);
  fprintf(fout, "{nullptr, nullptr, 0, nullptr}};\n\n");

  fprintf(fout, "#ifdef VTK_PY3K\n");
  fprintf(fout, "static PyModuleDef Py%s_Module = {\n", libName);
  fprintf(fout, "  PyModuleDef_HEAD_INIT,\n");
  fprintf(fout, "  \"%s\", // m_name\n", libName);
  fprintf(fout, "  nullptr, // m_doc\n");
  fprintf(fout, "  0, // m_size\n");
  fprintf(fout, "  Py%s_Methods, //m_methods\n", libName);
  fprintf(fout, "  nullptr, // m_reload\n");
  fprintf(fout, "  nullptr, // m_traverse\n");
  fprintf(fout, "  nullptr, // m_clear\n");
  fprintf(fout, "  nullptr  // m_free\n");
  fprintf(fout, "};\n");
  fprintf(fout, "#endif\n\n");

  fprintf(fout, "extern  \"C\" {%sPyObject *real_init%s(const char * /*unused*/); }\n\n", dllexp,
    libName);

  fprintf(fout, "PyObject *real_init%s(const char * /*unused*/)\n{\n", libName);

  /* module init function */
  fprintf(fout, "#ifdef VTK_PY3K\n");
  fprintf(fout, "  PyObject *m = PyModule_Create(&Py%s_Module);\n", libName);
  fprintf(fout, "#else\n");
  fprintf(fout,
    "  PyObject *m = Py_InitModule(\"%s\",\n"
    "                              Py%s_Methods);\n",
    importName, libName);
  fprintf(fout, "#endif\n\n");

  fprintf(fout, "  PyObject *d = PyModule_GetDict(m);\n");
  fprintf(fout, "  if (!d)\n");
  fprintf(fout, "  {\n");
  fprintf(fout, "    Py_FatalError(\"can't get dictionary for module %s\");\n", libName);
  fprintf(fout, "  }\n\n");

  /* import all the modules that we depend on */
  if (numDepends > 0)
  {
    fprintf(fout, "  const char *depends[%d] = {\n", numDepends);
    for (i = 0; i < numDepends; i++)
    {
      fprintf(fout, "    \"%s\",\n", depends[i]);
    }
    fprintf(fout, "  };\n\n");
    fprintf(fout, "  for (int i = 0; i < %d; i++)\n", numDepends);
    fprintf(fout, "  {\n");
    fprintf(fout, "    if (!vtkPythonUtil::ImportModule(depends[i], d))\n");
    fprintf(fout, "    {\n");
    fprintf(fout,
      "      return PyErr_Format(PyExc_ImportError,\n"
      "        \"Failed to load %s: No module named %%s\",\n"
      "        depends[i]);\n",
      libName);
    fprintf(fout, "    }\n");
    fprintf(fout, "  }\n\n");
  }

  for (i = 0; i < numFiles; i++)
  {
    fprintf(fout, "  PyVTKAddFile_%s(d);\n", files[i]);
  }

  fprintf(fout, "\n");
  fprintf(fout, "  vtkPythonUtil::AddModule(\"%s\");\n\n", libName);

  fprintf(fout, "  return m;\n");
  fprintf(fout, "}\n\n");
}

int main(int argc, char* argv[])
{
  FILE* file;
  FILE* fout_init;
  FILE* fout_impl;
  int numFiles = 0;
  int numDepends = 0;
  char libName[250];
  char importName[250];
  char tmpVal[250];
  char* files[4000];
  char* depends[400];
  int doDepends = 0;

  if (argc < 4)
  {
    fprintf(stderr, "Usage: %s input_file init_file impl_file [optional prefix]\n", argv[0]);
    return 1;
  }

  file = fopen(argv[1], "r");
  if (!file)
  {
    fprintf(stderr, "Input file %s could not be opened\n", argv[1]);
    return 1;
  }

  /* read the info from the file */
  if (fscanf(file, "%249s", libName) != 1)
  {
    fprintf(stderr, "Error getting libName\n");
    fclose(file);
    return 1;
  }
  /* read in the classes */
  while (fscanf(file, "%249s", tmpVal) != EOF)
  {
    if (strcmp(tmpVal, "DEPENDS") == 0)
    {
      doDepends = 1;
    }
    else if (doDepends)
    {
      depends[numDepends++] = strdup(tmpVal);
    }
    else
    {
      files[numFiles++] = strdup(tmpVal);
    }
  }
  /* close the file */
  fclose(file);
  file = NULL;

  fout_init = fopen(argv[2], "w");
  if (!fout_init)
  {
    return 1;
  }

  fout_impl = fopen(argv[3], "w");
  if (!fout_impl)
  {
    fclose(fout_init);
    return 1;
  }

  if (argc == 5)
  {
    size_t prefix_len = strlen(argv[4]);
    size_t lib_len = strlen(libName);
    memcpy(importName, argv[4], prefix_len);
    memcpy(importName + prefix_len, libName, lib_len);
    importName[prefix_len + lib_len] = '\0';
  }
  else
  {
    strcpy(importName, libName);
  }

  /* extra functions, types, etc. for the CommonCore module */
  if (strcmp(libName, "vtkCommonCore") == 0 || strcmp(libName, "vtkCommonKit") == 0)
  {
    files[numFiles] = strdup("PyVTKExtras");
    numFiles++;
  }

  CreateInitFile(libName, fout_init);
  CreateImplFile(libName, importName, numDepends, depends, numFiles, files, fout_impl);
  fclose(fout_init);
  fclose(fout_impl);

  return 0;
}