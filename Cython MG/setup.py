#!/usr/bin/env python
"""
Created on Thu Mar  2 23:49:04 2017

@author: Admin
"""

from distutils.core import setup
from Cython.Build import cythonize

setup(
  ext_modules = cythonize(["Verlet_main_MG.pyx"]),
)
