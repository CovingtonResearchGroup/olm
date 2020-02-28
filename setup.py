#!/usr/bin/env python
from distutils.core import setup
from Cython.Build import cythonize

setup(name='olm',
    version='0.35',
    author='Matt Covington',
    url='https://github.com/speleophysics/olm/',
    packages=['olm', 'olm.USGS', 'olm.loggers'],
    ext_modules = cythonize('olm/calcite.pyx')
    )
    
