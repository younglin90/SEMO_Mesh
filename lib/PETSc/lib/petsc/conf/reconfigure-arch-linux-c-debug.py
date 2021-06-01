#!/usr/bin/python3
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    '--with-cc=mpiicc',
    '--with-cxx=mpiicpc',
    '--with-fc=mpiifort',
    'PETSC_ARCH=arch-linux-c-debug',
  ]
  configure.petsc_configure(configure_options)
