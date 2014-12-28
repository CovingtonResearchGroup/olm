waterchem.USGS package
======================

This package contains functions for retrieving and processing data from the USGS National Water Information System database (NWIS).  Using these modules it is possible to retrieve water quality data, discharge data, and site metadata.  These data are then converted into Pandas data objects and stored in pickled and spreadsheet format for later usage. The package also contains functions for processing the Pandas data objects through the USGS water chemistry code PHREEQC.  For this, an installation of PHREEQC is required.

The most typical entry point is to use the convenience function ``runWQXtoPandas``.  This function can be used from within python by importing it from the ``waterchem.USGS.WQXtoPandas`` module.  It can additionally be used from the command line using::
    
    $./WQXtoPandas.py <name of input xls file>

Most functionality can be controlled from within the Excel file.  A template is provided with the waterchem package in the examples directory. The examples file also contains some brief code illustrating how data can be read in from the pickled outputs of ``WQXtoPandas`` and plotted.


:mod:`waterchem.USGS.WQXtoPandas` module
----------------------------------------

.. currentmodule:: waterchem.USGS.WQXtoPandas

.. autosummary::
   :toctree:

   runWQXtoPandas
   WQXtoPandas

:mod:`waterchem.USGS.PhreeqcPandas` module
------------------------------------------

.. currentmodule:: waterchem.USGS.PhreeqcPandas

.. autosummary::
   :toctree:

   processPanel
   readPhreeqcOutput
   writePhreeqcInput

:mod:`waterchem.USGS.DataRetrieval` module
------------------------------------------

.. currentmodule:: waterchem.USGS.DataRetrieval

.. autosummary::
   :toctree:

   GetDailyDischarge
   GetSiteData
   querySiteList

:mod:`waterchem.USGS.siteListExtraction` module
-----------------------------------------------

.. currentmodule:: waterchem.USGS.siteListExtraction

.. autosummary::
   :toctree:

   extractSitesFromText
   extractSitesFromXML	


