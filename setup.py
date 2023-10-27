#!/usr/bin/env python
from distutils.core import setup


setup(name='olm-karst',
    version='0.40',
    author='Matthew D. Covington',
      author_email = 'mcoving@uark.edu',
    url='https://github.com/speleophysics/olm/',
      download_url = 'https://github.com/CovingtonResearchGroup/olm/archive/refs/tags/0.40.tar.gz',
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
    
