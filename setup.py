#!/usr/bin/env python
from distutils.core import setup


setup(name='olm',
    version='0.37',
    author='Matthew D. Covington',
      author_email = 'mcoving@uark.edu',
    url='https://github.com/speleophysics/olm/',
      download_url = 'https://github.com/CovingtonResearchGroup/olm/archive/refs/tags/0.37.tar.gz',
    packages=['olm', 'olm.USGS', 'olm.loggers'],
      install_requires = ['numpy',
                          'scipy',
                          'matplotlib',
                          'pandas',
                          'lxml',
                          'xlrd',
                          'xlwt',
                          'openpyxl',
                          'requests',
      ],
      
    )
    
