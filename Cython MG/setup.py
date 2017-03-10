#!/usr/bin/env python
"""
Created on Thu Mar  2 23:49:04 2017

@author: Admin
"""
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

"""
setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [
        Extension("Verlet_main_MG",
                  ["Verlet_main_MG.pyx"],
                  include_dirs=[numpy.get_includes()],
                  libraries=["m"],
        ),
])"""

setup(
  ext_modules = cythonize(["Verlet_IC_MG.pyx"]),
)
