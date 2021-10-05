olm
===

A python package for common water chemistry calculations and automated processing of water quality data. The package is specifically oriented toward common carbonate equilibrium chemistry calculations and calculation of calcite dissolution rates. It is possible to process data from a variety of sources, including time series data from data loggers. However, a variety of functions have already been incorporated that allow automated retrieval of USGS water quality data that meet specified criteria. These data can also be automatically run through speciation calculations using PHREEQC. Additional analysis and plotting functions are currently being developed.

Installation
------------
**Installing using pip**
```sh
pip install olm-karst
```


**Installing from source**

Unzip the source code and run:
```sh
python setup.py install
```
However, you will first want to make sure you have the following dependencies installed.

Dependencies
------------
 * NumPy
 * SciPy
 * matplotlib
 * pandas
 * lxml
 * xlrd
 * xlwt
 * openpyxl
 * requests

To compile the documentation, you will also need 'numpydoc'. This package should now be compatible with Python 2 or 3. If you have trouble in either of these Python versions, please create an issue that explains the problem.

Full documentation can be found at:

http://olm.readthedocs.org/

and up to date versions of the code are available at:

http://www.github.com/speleophysics/olm

This package is a work in progress and undergoes frequent updates. If you find that something is not working as expected, a potential bug, or have suggestions for features that would be useful to you, feel free to submit such comments as issues on Github. Alternatively, you can also find my contact information on my webpage:

http://www.speleophysics.com
