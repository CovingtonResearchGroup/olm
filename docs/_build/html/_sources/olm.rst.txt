olm Package
=================

The olm package provides a chemical ``solution`` object along with a suite of functions for conducting geochemical calculations.  The package is currently predominantly oriented toward calcium carbonate geochemistry and common calculations conducted for karst waters.  The ``calcite`` module includes functions for creating solutions from a few chemical parameters (such as pCO2 and hardness), as well as functions for calculating equilibrium solutions and dissolution rates according to the PWP equation. Documentation and code development are still ongoing.

:mod:`olm.calcite` Module
-------------------------------

.. currentmodule:: olm.calcite

.. autosummary::
   :toctree:
   
   solutionFromCaPCO2
   solutionFromCaPCO2Relaxed
   solutionFrompHCaRelaxed
   H2CO3fromPCO2
   H2CO3sfromPCO2
   PCO2FromSolution
   activityHFromPCO2
   concCaEqFromPCO2
   PCO2EqFromCa
   concCaEqFromSolution
   concHFromCaPCO2Relaxed
   pwpFromSolution
   pwpRateTheory
   pwpRatePascal
   pwpRateFranci
   palmerRate
   palmerFromSolution
   createPalmerInterpolationFunctions
   dissRateFromCaPCO2
   pwp_to_mm_yr
   calc_K_H
   calc_K_W
   calc_K_0
   calc_K_1
   calc_K_2
   calc_K_c
   calc_kappa1
   calc_kappa2
   calc_kappa3
   calc_kappa4Theory	
   calc_kappa4Pascal
   calc_kappa4Franci
   calc_k1
   calc_k2
   calc_k_neg1
   calc_k_neg2


:mod:`olm.general` Module
-------------------------------

.. currentmodule:: olm.general

.. autosummary:: 
   :toctree: 

   solution
   condTo25
   HardnessFromCond
   CaFromCond
   CtoK
   KtoC
   DebyeHuckel
   neutralGamma
   gamma
   approxI
   molL_to_mgL
   mgL_to_molL
   mmolL_to_meqL
   molL_to_meqL
   getProperties

Subpackages
-----------

.. toctree::
   
   olm.USGS
   olm.loggers
