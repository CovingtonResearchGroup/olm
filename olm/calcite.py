#####################################################################
## Functions related to water chemistry and dissolution of Calcite ##
#####################################################################

import numpy as np
import pandas
from general import *
from scipy.optimize import brentq#, fminbound
from scipy.optimize import fsolve

#Define some useful constants
R = 8.3145 #J / (mol * K)
#H2OmolPerL = 55.5

def PCO2FromSolution(sol):
    """
    Calculate partial pressure of CO2 from a solution object.

    Parameters
    ----------
    sol : solution object, numpy.ndarray of solution objects, or pandas Series of solution objects

    Returns
    -------
    pCO2 : float, numpy.ndarray, or pandas series
       partial pressure(s) of CO2 for the solution(s)

    Notes
    -----
    Assumes a H20-CO2-CaCO3 system.  Uses equation 2.30 from Dreybrodt (1988) and assumes an open system.

    """
    def calc_PCO2(this_sol):
        Ca_conc = this_sol.ions['Ca']['conc_mol']    
        gamma_H = this_sol.activity_coef('H')
        gamma_HCO3 = this_sol.activity_coef('HCO3')
        gamma_Ca = this_sol.activity_coef('Ca')
        pH = this_sol.pH
        H_conc = 10.0**(-pH)/gamma_H
        #calculate mass action constans 
        T_K = this_sol.T_K    
        K_c = calc_K_c(T_K)
        K_2 = calc_K_2(T_K)
        K_1 = calc_K_1(T_K)
        K_0 = calc_K_0(T_K)
        K_H = calc_K_H(T_K)
        #pCO2 derived from equation 2.30 in Dreybrodt 1988 and assuming an
        #open system, where f approaches infty.  See notebook for details.
        pCO2 = (gamma_H * gamma_HCO3 / (K_1*K_H*(1+1/K_0)) ) * \
               (H_conc**2. + 2.0*H_conc*Ca_conc)
        return pCO2
    is_series = (type(sol)==pandas.core.series.Series) or (type(sol)==pandas.core.series.TimeSeries)
    if (type(sol)==np.ndarray) or is_series:
        pCO2 = np.zeros(np.size(sol))
        for i, single_sol in enumerate(sol):
            pCO2[i] = calc_PCO2(single_sol)
        if is_series:
            pCO2 = pandas.Series(pCO2,index=sol.index)
    else:
        pCO2 = calc_PCO2(sol)
    return pCO2

def concCaEqFromSolution(sol):
    """
    Calculates the equilibrium concentration of calcium for a solution object.

    First calculates the partial pressure of CO2, and then uses PCO2 and temperature to calculate equilibrium Ca.

    Parameters
    ----------
    sol : solution object, numpy.ndarray of solution objects, or pandas Series of solution objects

    Returns
    -------
    CaEq : float, numpy.ndarray, or pandas Series
       Equilibrium concentration(s) of CaEq for the given solution(s) in mol/L.

    Notes
    -----
    Assumes a H20-CO2-CaCO3 system.
    """
    PCO2 = PCO2FromSolution(sol)
    is_series = (type(sol)==pandas.core.series.Series) or (type(sol)==pandas.core.series.TimeSeries)
    if (type(sol)==np.ndarray) or is_series:
        T_C = np.zeros(np.size(sol))
        for i, single_sol in enumerate(sol):
            T_C[i] = single_sol.T_C
    else:
        T_C = sol.T_C
    CaEq = concCaEqFromPCO2(PCO2, T_C = T_C)
    return CaEq

#Calculate the equilibrium concentration of Ca using PCO2 and T_C
def concCaEqFromPCO2(PCO2, T_C = 25.):
    """
    Calculates the equilibrium concentration of calcium using PCO2 and temp.

    Iteratively solves for the equilibrium concentration of calcium from PCO2 and temp.  First guesses that activity coefficients are 1, and then iterates to solution using scipy's brentq function. 

    Parameters
    ----------
    PCO2 : float, numpy.ndarray, or pandas Series
       partial pressure of CO2 (atm)
    T_C : float, numpy.ndarray, or pandas Series (optional)
       temperature of solution in degrees Celsius (default = 25 C)

    Returns
    -------
    CaEq : float, numpy.ndarray, or pandas Series
       Equilibrium concentration(s) of calcium in mol/L
 
    Notes
    -----
    Assumes a H20-CO2-CaCO3 system.   
    If a numpy array or pandas Series object are passed in as PCO2 arguments, then equilibrium concentrations will be found iteratively in a for-loop and returned as the same data type given in the argument.
    """
    
    def Ca_minimize(Ca,T_C_func,K_c,K_2,K_1,K_H,PCO2_func):
        if Ca<0:
            return 10.
        I = approxI(Ca)
        properties = getProperties()
        z_Ca = properties['Ca']['charge']
        r_Ca = properties['Ca']['radius']
        gamma_Ca = DebyeHuckel(I,z_Ca,r_Ca,T=T_C_func)
        z_HCO3 = properties['HCO3']['charge']
        r_HCO3 = properties['HCO3']['radius']
        gamma_HCO3 = DebyeHuckel(I,z_HCO3,r_HCO3,T=T_C_func)
        return Ca - (PCO2_func*K_1*K_c*K_H/(4.*K_2*gamma_Ca*gamma_HCO3**2.))**(1./3.)
    #If there is only one value for T_C, but multiple PCO2 values, make an array of equal values
    if (np.size(PCO2)>1) and (np.size(T_C) == 1):
        T_C = T_C + np.zeros(np.size(PCO2))
    T_K = CtoK(T_C)
    K_c = calc_K_c(T_K)
    K_2 = calc_K_2(T_K)
    K_1 = calc_K_1(T_K)
    K_H = calc_K_H(T_K)
    #make a guess assuming activities are = 1
    guess_Ca = (PCO2*K_1*K_c*K_H/(4.*K_2))**(1./3.)
    maxCa = 10.*guess_Ca
    minCa = guess_Ca
    is_series = (type(PCO2)==pandas.core.series.Series) or (type(PCO2)==pandas.core.series.TimeSeries)
    if (type(PCO2)==np.ndarray) or is_series:
        #We have a numpy array or pandas Series. Loop through solutions.
        CaEq = np.zeros(np.size(PCO2))       
        for i, single_PCO2 in enumerate(PCO2):
            try:
                CaEq[i] = brentq(Ca_minimize, guess_Ca[i], 10.*guess_Ca[i], args=(T_C[i],K_c[i],K_2[i],K_1[i],K_H[i], PCO2[i]))
            except RuntimeError:
                CaEq[i] = np.nan
        if is_series:
            #Create a pandas series object from the CaEq array
            CaEq = pandas.Series(CaEq, index=PCO2.index)
    else: #We only have a single value
        try:
            CaEq = brentq(Ca_minimize, guess_Ca, 10.*guess_Ca, args=(T_C,K_c,K_2,K_1,K_H,PCO2))
        except RuntimeError:
            CaEq = np.nan
    return CaEq

#Calculates equilibrium activity of H+ given PCO2
#  - uses relaxed charge balance assumption
def activityHFromPCO2(PCO2, T_C = 25., CaEq = None):
    """
    Calculates equilibrium activity of H+ given PCO2 using relaxed charge balance.

    Calculates hydrogen activity at equilibrium given PCO2, temperature, and (optionally) equilibrium calcium concentration (mol/L).  Assumes a relaxed charge balance (see 2.18a in Dreybrodt [1988]).  If keyword CaEq is not given, then it is iteratively calculated using concCaEqFromPCO2().

    Parameters
    ----------
    PCO2 : float
       partial pressure of CO2 (atm)
    T_C : float, optional
       temperature of solution in degrees Celsius (default = 25 C)
    CaEq : float
       Equilibrium calcium concentration (mol/L), optional

    Returns
    -------
    aHeq : float
       equilibrium activity of hydrogen ion (mol/L)

    Notes
    -----
    Assumes a H20-CO2-CaCO3 system.
    """
    if CaEq == None:
        CaEq = concCaEqFromPCO2(PCO2, T_C=T_C)
    I = approxI(CaEq)
    properties = getProperties()
    z_Ca = properties['Ca']['charge']
    r_Ca = properties['Ca']['radius']
    gamma_Ca = DebyeHuckel(I,z_Ca,r_Ca,T=T_C)
    z_HCO3 = properties['HCO3']['charge']
    r_HCO3 = properties['HCO3']['radius']
    gamma_HCO3 = DebyeHuckel(I,z_HCO3,r_HCO3,T=T_C)
    T_K = CtoK(T_C)
    K_c = calc_K_c(T_K)
    K_2 = calc_K_2(T_K)
    K_1 = calc_K_1(T_K)
    K_H = calc_K_H(T_K)
    a_Heq = ( ((K_1*K_H*PCO2)**2.)*K_2*gamma_Ca / (2*K_c*gamma_HCO3) )**(1./3.)
    return a_Heq

#Calculate H+ concentration from Ca and PCO2 assuming relaxed charge balance
def concHFromCaPCO2Relaxed(Ca, PCO2, T_C = 25.):
    """
    Calculates concentration of H+ from calcium concentration and PCO2 using relaxed charge balance.  Uses equation 2.30a from Dreybrodt (1988).

    Parameters
    ----------
    Ca : float
       concentration of calcium in mol/L
    PCO2 : float
       partial pressure of CO2 (atm)
    T_C : float, optional
       temperature of solution in degrees Celsius (default = 25 C)

    Returns
    -------
    concH : float
       concentration of hydrogen ions

    Notes
    -----
    Assumes a H20-CO2-CaCO3 system.
    """
    #from eqn 2.30a in Dreybrodt 1988
    T_K = CtoK(T_C)
    I = 3.*Ca
    gamma_H = gamma('H', I, T_C=T_C)
    gamma_HCO3 = gamma('HCO3', I, T_C=T_C)
    K_1 = calc_K_1(T_K)
    K_H = calc_K_H(T_K)
    K_0 = calc_K_0(T_K)
    HCO3_sqaured = K_1*K_H*PCO2*(1.+1./K_0)/(gamma_H*gamma_HCO3)
    concH = -Ca + np.sqrt(Ca**2. + HCO3_sqaured)
    return concH

#Calculate H+ concentration given Ca and PCO2, makes no relaxed charge
#balance assumption
def solutionFromCaPCO2(Ca, PCO2, T_C = 25., per_tol = 0.001):
    """
    Creates a solution object from a given concentration of calcium and PCO2.

    Parameters
    ----------
    Ca : float, numpy.ndarray, or pandas Series 
       concentration of calcium in mol/L
    PCO2 : float, numpy.ndarray, or pandas Series 
       partial pressure of CO2 (atm)
    T_C : float, , numpy.ndarray, or pandas Series (optional)
       temperature of solution in degrees Celsius (default = 25 C)
    per_tol : float
       the fractional change in H concentration between iterations upon which the iteration is terminated
    

    Returns
    -------
    sol : solution object, numpy.ndarray of solution objects, or pandas Series of solution objects       

    Notes
    -----
    Assumes a H20-CO2-CaCO3 system.  Guesses concentration of H using relaxed charge balance assumption, and then iterates to full solution.
    """
    def calc_sol(Ca_in,PCO2_in,T_C_in,per_tol=0.001):
        I_guess = 3.*Ca_in
        T_K = CtoK(T_C_in)
        H_guess = concHFromCaPCO2Relaxed(Ca_in,PCO2_in, T_C=T_C_in)
        K_W = calc_K_W(T_K)
        K_H = calc_K_H(T_K)
        K_0 = calc_K_0(T_K)
        K_1 = calc_K_1(T_K)
        K_2 = calc_K_2(T_K)
        K_6 = K_1*(1.+1./K_0)
        found=False
        while not(found):
            #estimate activity coefficients
            gamma_H = gamma('H', I_guess, T_C=T_C_in)
            gamma_OH = gamma('OH', I_guess, T_C=T_C_in)
            gamma_HCO3 = gamma('HCO3', I_guess, T_C=T_C_in)
            gamma_CO3 = gamma('CO3', I_guess, T_C=T_C_in)        
            H_new= fsolve(lambda H: 2.*Ca_in + H - K_W/(gamma_H*gamma_OH*H)\
                          - K_6*K_H*PCO2_in/(gamma_HCO3*gamma_H*H)*\
                          (1. + 2.*K_2*gamma_HCO3/(gamma_CO3*gamma_H*H)),\
                          H_guess)[0]
            #calculate ion concentrations from guess H+ concentration
            OH = K_W/(gamma_OH*gamma_H*H_new)
            CO3 = K_2*K_6*K_H*PCO2_in/((gamma_H*H_new)**2)/gamma_CO3
            HCO3 = K_6*K_H*PCO2_in/(gamma_HCO3*gamma_H*H_new)
            I_new = 0.5*(H_new + OH + HCO3 + 4.*CO3 + 4.*Ca_in)
            if (np.abs(H_new - H_guess)/H_guess < per_tol):
                found = True
            else:
                H_guess = H_new
                I_guess = I_new
            #calculate non-charge ions
        CO2 = K_H*PCO2_in
        H2CO3 = (K_H/K_0)*PCO2_in
        H2CO3s = H2CO3 + CO2
        pH = -np.log10(H_new*gamma_H)
        #creation solution with these species
        sol = solution(['H', 'OH', 'CO3', 'HCO3', 'Ca', 'CO2', 'H2CO3', 'H2CO3s'],
                       [H_new, OH, CO3, HCO3, Ca_in, CO2, H2CO3, H2CO3s],
                       "mol/L", T=T_C_in, pH = pH)
        return sol
    is_series = (type(Ca)==pandas.core.series.Series) or (type(Ca)==pandas.core.series.TimeSeries)
    if (type(Ca)==np.ndarray) or is_series:
        sol_arr = np.empty(np.size(Ca),dtype=object)        
        for i, this_Ca in enumerate(Ca):
            if np.size(T_C)==1:
                sol_arr[i] = calc_sol(Ca[i],PCO2[i],T_C,per_tol=per_tol)
            else:
                sol_arr[i] = calc_sol(Ca[i],PCO2[i],T_C[i],per_tol=per_tol)

        if is_series:
            sol_arr = pandas.Series(sol_arr,index=Ca.index)
        return sol_arr
    else:
        return calc_sol(Ca,PCO2,T_C,per_tol=per_tol)

        
# Function to calculate H+ concentration from Calcium concentration and pC02 
# using approximation in equation 2.30a (with an additional assumption that 
# chi -> infty, as for an open system) from Dreybrodt 1988. which assumes
# pH < 8 such that carbonate and OH- species can be neglected
#   Ca = Calcium concentration mol/L
#   PCO2 = partial pressure of CO2
def solutionFromCaPCO2Relaxed(Ca, PCO2, T_C = 25.):
    """
    Creates a solution object from a given concentration of calcium and PCO2.

    Creates a solution object from a given concentration of calcium, PCO2, and optional temperature.  Uses the approximate charge balance assumption (equation 2.30a in Dreybrodt [1988]).  This is valid when pH < 8, such that CO3- and OH- species can be neglected.  

    Parameters
    ----------
    Ca : float
       concentration of calcium in mol/L
    PCO2 : float
       partial pressure of CO2 (atm)
    T_C : float, optional
       temperature of solution in degrees Celsius (default = 25 C)

    Returns
    -------
    sol : solution object
       
    Notes
    -----
    Assumes a H20-CO2-CaCO3 system.  Calculates concentration of H using relaxed charge balance assumption.
    """
    T_K = CtoK(T_C)
    H2CO3s = H2CO3sfromPCO2(PCO2, T_K=T_K)
    H2CO3 = H2CO3fromPCO2(PCO2, T_K=T_K)
    properties = getProperties()
    I = approxI(Ca)
    gamma_H = DebyeHuckel(I, 
                          properties['H']['charge'],
                          properties['H']['radius'],
                          T = T_C)
    gamma_HCO3 = DebyeHuckel(I, 
                          properties['HCO3']['charge'],
                          properties['HCO3']['radius'],
                          T = T_C)
    gamma_CO3 = DebyeHuckel(I,
                            properties['CO3']['charge'],
                            properties['CO3']['radius'],
                            T = T_C)
    K_0 = calc_K_0(T_K)
    K_1 = calc_K_1(T_K)
    K_2 = calc_K_2(T_K)
    K_H = calc_K_H(T_K)
    K_6 = K_1*(1.+1./K_0)
    v_over_w = K_1*H2CO3s/(gamma_H*gamma_HCO3)
    #equation 2.30a with chi-->infty
    H = -Ca + 0.5*np.sqrt(4.*(Ca**2.)+ 4.*v_over_w)
    #from equation 2.24 with OH and CO3 neglected
    HCO3 = 2.*Ca + H
    CO3 = K_2*K_6*K_H*PCO2/((gamma_H*H)**2)/gamma_CO3
    #calculate amount fraction of H (neglecting other ions)
    H_activity = gamma_H*H
    pH = -np.log10(H_activity)
    sol = solution(['H', 'Ca', 'HCO3', 'H2CO3', 'H2CO3s', 'CO3'], [H, Ca, HCO3, H2CO3, H2CO3s, CO3], units="mol/L", T=T_C, T_units='C', pH=pH)
    return sol

# Function to calculate solution from Calcium concentration and pH
# using approximation in equation 2.30a (with an additional assumption that 
# chi -> infty, as for an open system) from Dreybrodt 1988. which assumes
# pH < 8 such that carbonate and OH- species can be neglected
#   Ca = Calcium concentration mol/L
#   pH = partial pressure of CO2
def solutionFrompHCaRelaxed(Ca, pH, T_C = 25.):
    """
    Creates a solution object from a given concentration of calcium and pH.

    Creates a solution object from a given concentration of calcium, pH, and optional temperature.  Uses the approximate charge balance assumption (equation 2.30a in Dreybrodt [1988]).  This is valid when pH < 8, such that CO3- and OH- species can be neglected.  

    Parameters
    ----------
    Ca : float
       concentration of calcium in mol/L
    pH : float
       pH
    T_C : float, optional
       temperature of solution in degrees Celsius (default = 25 C)

    Returns
    -------
    sol : solution object
       
    Notes
    -----
    Assumes a H20-CO2-CaCO3 system.  Calculates concentration of H using relaxed charge balance assumption.
    """
    T_K = CtoK(T_C)
    properties = getProperties()
    I = approxI(Ca)
    gamma_H = DebyeHuckel(I, 
                          properties['H']['charge'],
                          properties['H']['radius'],
                          T = T_C)
    gamma_HCO3 = DebyeHuckel(I, 
                          properties['HCO3']['charge'],
                          properties['HCO3']['radius'],
                          T = T_C)
    gamma_CO3 = DebyeHuckel(I,
                            properties['CO3']['charge'],
                            properties['CO3']['radius'],
                            T = T_C)
    
    gamma_Ca = DebyeHuckel(I, 
                           properties['Ca']['charge'],
                           properties['Ca']['radius'],
                           T = T_C)
    H = 10.0**(-pH)/gamma_H
    #calculate mass action constans 
    K_c = calc_K_c(T_K)
    K_2 = calc_K_2(T_K)
    K_1 = calc_K_1(T_K)
    K_0 = calc_K_0(T_K)
    K_H = calc_K_H(T_K)
    K_6 = K_1*(1.+1./K_0)
    #pCO2 derived from equation 2.30 in Dreybrodt 1988 and assuming an
    #open system, where f approaches infty.  See notebook for details.
    pCO2 = (gamma_H * gamma_HCO3 / (K_1*K_H*(1+1/K_0)) ) * \
        (H**2. + 2.0*H*Ca)
    H2CO3 = H2CO3fromPCO2(pCO2, T_K=T_K)
    H2CO3s = H2CO3sfromPCO2(pCO2, T_K=T_K)
    CO3 = K_2*K_6*K_H*pCO2/((gamma_H*H)**2)/gamma_CO3
    K_1 = calc_K_1(T_K)
    #equation 2.30a with chi-->infty
    HCO3 = 2.*Ca + H
    is_series = (type(Ca)==pandas.core.series.Series) or (type(Ca)==pandas.core.series.TimeSeries)
    if (type(Ca)==np.ndarray) or is_series:
        sol_arr = np.empty(np.size(Ca),dtype=object)
        for i in range(np.size(Ca)):
            if np.size(T_C)==1:
                sol_arr[i] = solution(['H', 'Ca', 'HCO3', 'H2CO3', 'H2CO3s', 'CO3'], [H[i], Ca[i], HCO3[i], H2CO3[i], H2CO3s[i], CO3[i]], units="mol/L", T=T_C, T_units='C', pH=pH[i])
            else:
                sol_arr[i] = solution(['H', 'Ca', 'HCO3', 'H2CO3', 'H2CO3s', 'CO3'], [H[i], Ca[i], HCO3[i], H2CO3[i], H2CO3s[i], CO3[i]], units="mol/L", T=T_C[i], T_units='C', pH=pH[i])
        if is_series:
            sol_arr = pandas.Series(sol_arr, index=Ca.index)
        return sol_arr
    else:
        sol = solution(['H', 'Ca', 'HCO3', 'H2CO3', 'H2CO3s', 'CO3'], [H, Ca, HCO3, H2CO3, H2CO3s, CO3], units="mol/L", T=T_C, T_units='C', pH=pH)
        return sol


#Calculate concentration of Carbonic acid in equilibrium with a certain pCO2
def H2CO3fromPCO2(PCO2, T_K = 273.15 + 25., T_C = None):
    """
    Calculate concentration of carbonic acid in equilibrium with a certain PCO2.
    
    Parameters
    ----------
    PCO2 : float
       partial pressure of CO2
    T_K : float, optional       temperature in degrees Kelvin (default = 273.15 + 25)
    T_C : float, optional
       temperature in degrees Celsius (default = None).  If None, function assumes T_K was given or uses default value.  If T_C is given in function call, then function uses T_C value to calculate T_K.
       
    Returns
    -------
    H2CO3 : float
       concentration of carbonic acid in mol/L
    """
    if T_C != None:
        T_K = CtoK(T_C)
    K_H = calc_K_H(T_K)
    K_0 = calc_K_0(T_K)    
    H2CO3 = K_H*PCO2/K_0
    return H2CO3

#Calculate concentration of Carbonic acid in equilibrium with a certain pCO2
def H2CO3sfromPCO2(PCO2, T_K = 273.15 + 25., T_C = None):
    """
    Calculate concentration of carbonic acid + aqueous CO2 in equilibrium with a certain PCO2. [H2CO3s] = [H2CO3] + [CO2]
    
    Parameters
    ----------
    PCO2 : float
       partial pressure of CO2
    T_K : float, optional       temperature in degrees Kelvin (default = 273.15 + 25)
    T_C : float, optional
       temperature in degrees Celsius (default = None).  If None, function assumes T_K was given or uses default value.  If T_C is given in function call, then function uses T_C value to calculate T_K.
       
    Returns
    -------
    H2CO3s : float
       concentration of carbonic acid + aqueous CO2 in mol/L
    """
    if T_C != None:
        T_K = CtoK(T_C)
    K_H = calc_K_H(T_K)
    K_0 = calc_K_0(T_K)    
    H2CO3s = K_H*PCO2*(1+1/K_0)
    return H2CO3s

def pwpFromSolution(sol, PCO2=None, method='theory'):
    """
    Calculates the PWP dissolution rate from a solution object.

    Parameters
    ----------
    sol : solution object, numpy.ndarray, or pandas Series
       An olm solution object for which the calcite dissolution rate will be calculated.
    PCO2 : float
       The partial pressure of CO2 for the solution.  If not given, it will be calculated from the solution object using PCO2FromSolution().
    method : string
       A string that is equal to 'theory', 'pascal', or 'franci' that specifies the version of the PWP equation to use.

    Returns
    -------
    R : float, numpy.ndarray, or pandas Series
       calcite dissolution rate according to the PWP equation (mmol/cm^2/s)
       
    """
    if PCO2==None:
        PCO2 = PCO2FromSolution(sol)

    def calc_rate(sol_in,PCO2_in):    
        #Check whether all necessary ions are present
        if ('Ca' in sol_in.ions) and ('H' in sol_in.ions) and ('HCO3' in sol_in.ions):
            #Check whether H2CO3s is present, and calculate if necessary
            if not 'H2CO3s' in sol_in.ions:
                K_H = calc_K_H(sol_in.T_K)
                CO2 = K_H*PCO2_in
                if not 'H2CO3' in sol_in.ions:
                    K_0 = calc_K_0(sol_in.T_K)
                    a_H2CO3s = CO2*(1.+1./K_0)
                else:
                    a_H2CO3s = CO2 + sol_in.activity('H2CO3')
            else: #If we have it already, just read it
                a_H2CO3s = sol_in.activity('H2CO3s')
            #Pull out other ion concentrations
            a_Ca = sol_in.activity('Ca')
            a_H = sol_in.activity('H')
            a_HCO3 = sol_in.activity('HCO3')
            T_K = sol_in.T_K
            if method=='theory':
                R = pwpRateTheory(a_Ca=a_Ca, a_H2CO3s=a_H2CO3s, a_H=a_H, a_HCO3=a_HCO3, T_K=T_K,PCO2=PCO2_in)
            elif method=='pascal':
                R = pwpRatePascal(a_Ca=a_Ca, a_H2CO3s=a_H2CO3s, a_H=a_H, a_HCO3=a_HCO3, T_K=T_K,PCO2=PCO2_in)
            elif method=='franci':
                R = pwpRateFranci(a_Ca=a_Ca, a_H2CO3s=a_H2CO3s, a_H=a_H, a_HCO3=a_HCO3, T_K=T_K,PCO2=PCO2_in)
            else:
                print "method must be set to 'theory', 'pascal', or 'franci'"
                return -1
            return R
        else:
            print "Not all necessary ions are present in the solution object."
            return -1 
    is_series = (type(sol)==pandas.core.series.Series) or (type(sol)==pandas.core.series.TimeSeries)
    if (type(sol)==np.ndarray) or is_series:
        sol_arr = np.empty(np.size(sol),dtype=object)
        for i, this_sol in enumerate(sol):            
            sol_arr[i] = calc_rate(this_sol,PCO2[i])
        if is_series:
            sol_arr = pandas.Series(sol_arr, index=sol.index)
        return sol_arr
    else:
        return calc_rate(sol,PCO2)
        
#Calculate dissolution rate from PWP equation using an input concentrations
# kappa4 is calculated using relation from Dreybrodt's PASCAL code. 
def pwpRatePascal(a_Ca=0., a_H2CO3s=0., a_H=0., a_HCO3=0., T_K=25.+273.15,PCO2=0.):
    """
    Calculates PWP rate using relation for kappa4 found in PASCAL code.

    Calculates PWP dissolution rate for calcite using the relation for kappa4 that is found in the PASCAL code in Dreybrodt (1988).  This is also given in equation 3 from Buhmann and Dreybrodt (1985), The kinetics of calcite dissolution and precipitation in geologically relevant situations of karst areas: 1. Open system.  They say that it is a fit to the experimental data of Plummer et al. for values of PCO2 < 0.05 atm.

    Parameters
    ----------
    a_Ca : float
       activity of calcium (mol/L)
    a_H2CO3s : float
       activity of carbonic acid (mol/L)
    a_H : float
       activity of hydrogen (mol/L)
    a_HCO3 : float
       activity of bicarbonate (mol/L)
    T_K : float
       temperature degrees Kelvin
    PCO2 : float
       partial pressure of CO2 (atm)

    Returns
    -------
    R : float
       calcite dissolution rate according to the PWP equation (mmol/cm^2/s)

    """
    kappa1 = calc_kappa1(T_K)
    kappa2 = calc_kappa2(T_K)
    kappa3 = calc_kappa3(T_K)
    kappa4 = calc_kappa4Pascal(T_K,PCO2)
    R = kappa1*a_H + kappa2*a_H2CO3s + kappa3 - kappa4*a_Ca*a_HCO3
    return R

#Calculate dissolution rate from PWP equation using an input
# concentrations kappa4 is calculated using theoretical relation with
# a_H equal to current value, as done in Franci's code
def pwpRateFranci(a_Ca=0., a_H2CO3s=0., a_H=0., a_HCO3=0., T_K=25.+273.15,PCO2=0.):
    """
    Calculates PWP rate using relation for kappa4 used in Franci Gabrovsek's code.

    Calculates PWP rate using relation for kappa4 used in Franci Gabrovsek's speleogenesis code (pers. commun.).  This slight difference was discovered during testing of this code against Franci's calculations.  In this case, a_H in the equation for kappa4 is the bulk value and not the equilibrium surface value for the given carbonic acid concentration.

    Parameters
    ----------
    a_Ca : float
       activity of calcium (mol/L)
    a_H2CO3s : float
       activity of carbonic acid (mol/L)
    a_H : float
       activity of hydrogen (mol/L)
    a_HCO3 : float
       activity of bicarbonate (mol/L)
    T_K : float
       temperature degrees Kelvin
    PCO2 : float
       partial pressure of CO2 (atm)

    Returns
    -------
    R : float
       calcite dissolution rate according to the PWP equation (mmol/cm^2/s)
"""
    kappa1 = calc_kappa1(T_K)
    kappa2 = calc_kappa2(T_K)
    kappa3 = calc_kappa3(T_K)
    kappa4 = calc_kappa4Franci(T_K, a_H, a_H2CO3s)
    R = kappa1*a_H + kappa2*a_H2CO3s + kappa3 - kappa4*a_Ca*a_HCO3
    return R

def pwpRateTheory(a_Ca=0., a_H2CO3s=0., a_H=0., a_HCO3=0., T_K=25.+273.15,PCO2=0.):
    """
    Calculates PWP rate using theoretical relation for kappa4 from PWP.

    Calculates PWP rate using theoretical relation for kappa4 from PWP (as described by equation 25 in Plummer, Wigley, and Parkhurst (1978) and in Dreybrodt [1988] equation 6.22b).  In this case, a_H in the equation for kappa4 is the equilibrium surface value for the given carbonic acid concentration, as specified in the theory.

    Parameters
    ----------
    a_Ca : float
       activity of calcium (mol/L)
    a_H2CO3s : float
       activity of carbonic acid (mol/L)
    a_H : float
       activity of hydrogen (mol/L)
    a_HCO3 : float
       activity of bicarbonate (mol/L)
    T_K : float
       temperature degrees Kelvin
    PCO2 : float
       partial pressure of CO2 (atm)

    Returns
    -------
    R : float
       calcite dissolution rate according to the PWP equation (mmol/cm^2/s)
"""

    kappa1 = calc_kappa1(T_K)
    kappa2 = calc_kappa2(T_K)
    kappa3 = calc_kappa3(T_K)
    is_series = (type(a_Ca)==pandas.core.series.Series) or (type(a_Ca)==pandas.core.series.TimeSeries)
    if (type(a_Ca)==np.ndarray) or is_series:
        #We have a numpy array or pandas Series. Loop through and calculate rates individually    
        R = np.zeros(np.size(a_Ca))       
        for i, single_R in enumerate(R):
            kappa4 = calc_kappa4Theory(T_K[i], PCO2[i], a_H2CO3s[i])
            R[i] = kappa1[i]*a_H[i] + kappa2[i]*a_H2CO3s[i] + kappa3[i] - kappa4*a_Ca[i]*a_HCO3[i]
        if is_series:
            #Create a pandas series object from the R array
            R = pandas.Series(R, index=a_Ca.index)
    else:
        kappa4 = calc_kappa4Theory(T_K, PCO2, a_H2CO3s)
        R = kappa1*a_H + kappa2*a_H2CO3s + kappa3 - kappa4*a_Ca*a_HCO3
    return R

def pwp_to_mm_yr(R, rho=2.6):
    """
    Converts the PWP dissolution rates from mmol/cm^2/s to mm/year.

    Parameters
    ----------
    R : float
       Dissolution rate as provided by PWP rate functions in units of mmol/cm^2/s.
    rho : float
       Density of rock in g/cm^3.  Default is 2.6 g/cm^3, a typical value for limestone.

    Returns
    -------
    E : float
       Erosion rate in mm/year
    """
    #First convert from mmol to mol.
    R_mol = R*10.**(-3)
    properties = getProperties()
    CaCO3_weight = properties['Ca']['weight'] + properties['CO3']['weight']
    #Convert to grams
    R_g = R_mol*CaCO3_weight
    #R in cm/s
    R_cm_s = R_g/rho
    #Convert to mm
    R_mm_s = R_cm_s*10.
    #Convert from seconds to years
    E = R_mm_s*365.*24.*60.*60.
    return E

#Functions to calculate rate constants in PWP Equation from Equations
# 6.13, 6.14, and 6.22b in Dreybrodt 1988
def calc_kappa1(T_K):
    """
    Calculates kappa1 in the PWP equation.

    Calculates kappa1 in the PWP equation, according to Equation 5 in Plummer, Wigley, and Parkhurst (1978) or Equation 6.13 of Dreybrodt (1988).

    Parameters
    ----------
    T_K : float
       temperature Kelvin

    Returns
    -------
    kappa1 : float
       constant kappa1 in the PWP equation (cm/s)
    """
    kappa1 = 10.**(0.198 - 444./T_K)
    return kappa1

def calc_kappa2(T_K):
    """
    Calculates kappa2 in the PWP equation.

    Calculates kappa2 in the PWP equation, according to Equation 7 in Plummer, Wigley, and Parkhurst (1978) or Equation 6.14 of Dreybrodt (1988).

    Parameters
    ----------
    T_K : float
       temperature Kelvin

    Returns
    -------
    kappa2 : float
       constant kappa2 in the PWP equation (cm/s)
    """
    
    kappa2 = 10.**(2.84 - 2177./T_K)
    return kappa2
def calc_kappa3(T_K):
    """
    Calculates kappa3 in the PWP equation.

    Calculates kappa3 in the PWP equation, according to Equations 8 and 9 in Plummer, Wigley, and Parkhurst (1978) or Equations 6.14a and 6.14b of Dreybrodt (1988).

    Parameters
    ----------
    T_K : float
       temperature Kelvin

    Returns
    -------
    kappa3 : float
       constant kappa3 in the PWP equation (mmol/cm^2/s)
    """
    if (np.size(T_K)>1):
        kappa3 = np.zeros(np.size(T_K))
        for i, this_temp in enumerate(T_K):
            if (this_temp < 273.15+25):
                kappa3[i] = 10.**(-5.86 - 317./T_K[i])
            else:
                kappa3[i] = 10.**(-1.10 - 1737./T_K[i])

    else:
        if (T_K < 273.15+25):
            kappa3 = 10.**(-5.86 - 317./T_K)
        else:
            kappa3 = 10.**(-1.10 - 1737./T_K)
    return kappa3

def calc_kappa4Kaufmann(T_K, PCO2):
    """
    Calculates kappa4 in the PWP equation using the relation from Kaufmann and Dreybrodt (2007).

    Parameters
    ----------
    T_K : float
       temperature Kelvin
    PCO2 : float
       partial pressure of CO2 (atm)

    Returns
    -------
    kappa4 : float
       constant kappa4 in the PWP equation (cm^4/mmol/s)

    Notes
    -----
    
    """
    T_C = KtoC(T_K)
    if PCO2>0.05:
        kappa4 = 10.**(-2.375+0.025*T_C)
    else:
        kappa4 = 10.**(-2.375+0.025*T_C + 0.56*(-np.log10(PCO2)-1.3))
    return kappa4


def calc_kappa4Pascal(T_K,PCO2):
    """
    Calculates kappa4 in the PWP equation using fit from Buhmann and Dreybrodt (1985).

    Parameters
    ----------
    T_K : float
       temperature Kelvin
    PCO2 : float
       partial pressure of CO2 (atm)

    Returns
    -------
    kappa4 : float
       constant kappa4 in the PWP equation (cm^4/mmol/s)

    Notes
    -----
    See more info under documentation for pwpRatePascal().
    """
    T_C = KtoC(T_K)
    B = 3.077-0.0146*T_C
    kappa4 = 10.**(-B)*(1/PCO2)**0.611
    return kappa4
    
def calc_kappa4Franci(T_K, a_H, a_H2CO3s):
    """
    Calculates kappa4 in the PWP equation using approach from Franci's code.

    Parameters
    ----------
    T_K : float
       temperature Kelvin
    a_H : float
       activity of hydrogen (mol/L)
    a_H2CO3s : float
       activity of carbonic acid (mol/L)

    Returns
    -------
    kappa4 : float
       constant kappa4 in the PWP equation (cm^4/mmol/s)

    Notes
    -----
    See more info under documentation for pwpRateFranci().
    """
    K_2 = calc_K_2(T_K)
    K_c = calc_K_c(T_K)
    kappa1 = calc_kappa1(T_K)
    kappa2 = calc_kappa2(T_K)
    kappa3 = calc_kappa3(T_K)
    kappa4 = (K_2/K_c)*(kappa1 + 1/a_H*(kappa2*a_H2CO3s + kappa3) ) 
    return kappa4

def calc_kappa4Theory(T_K, PCO2, a_H2CO3s):
    """
    Calculates kappa4 in the PWP equation using the theoretical relation for kappa4 from Plummer, Wigley, and Parkhurst (1978) Equation 25 (as described in Dreybrodt [1988] equation 6.22b).  In this case, a_H in the equation for kappa4 is the equilibrium surface value for the given carbonic acid concentration, as specified in the theory.

    Parameters
    ----------
    T_K : float
       temperature Kelvin
    PCO2 : float
       partial pressure of CO2 (atm)
    a_H2CO3s : float
       activity of carbonic acid (mol/L)

    Returns
    -------
    kappa4 : float
       constant kappa4 in the PWP equation (cm/s)

    Notes
    -----
    See more info under documentation for pwpRateTheory().
    """
  
    T_C = KtoC(T_K)
    K_2 = calc_K_2(T_K)
    K_c = calc_K_c(T_K)
    kappa1 = calc_kappa1(T_K)
    kappa2 = calc_kappa2(T_K)
    kappa3 = calc_kappa3(T_K)
    #calculate equilbrium activity of H at surface
    a_Heq = activityHFromPCO2(PCO2, T_C=T_C)
    kappa4 = (K_2/K_c)*(kappa1 + 1/a_Heq*(kappa2*a_H2CO3s + kappa3) ) 
    return kappa4


#Functions for calculating mass action constants given temperature (K)
# from Dreybrodt, 1988
def calc_K_c(T_K):
    """
    Calculates equilibrium constant for calcite.

    Calculates equilibrium constant for calcite using equation from Table 2.2 in Dreybrodt (1988), originally reported in Plummer and Busenberg (1982).

    Parameters
    ----------
    T_K : float
       temperature Kelvin

    Returns
    -------
    K_c : float
       equilibrium constant for calcite
    """
    K_c = 10.**(-171.9065 - 0.077993*T_K + 2839.319/T_K + 71.595*np.log10(T_K))
    return K_c

def calc_K_2(T_K):
    """
    Calculates mass action constant for dissociation of bicarbonate.

    Calculates mass action constant for dissociation of bicarbonate using equation from Table 2.2 in Dreybrodt (1988), originally reported in Plummer and Busenberg (1982).

    Parameters
    ----------
    T_K : float
       temperature Kelvin

    Returns
    -------
    K_2 : float
       mass action constant for dissociation of bicarbonate
    """
    K_2 = 10.**(-107.8871 - 0.03252849*T_K + 5151.79/T_K + 38.92561*np.log10(T_K) - 563713.9/T_K/T_K)
    return K_2

def calc_K_1(T_K):
    """
    Calculates mass action constant for dissociation of carbonic acid.

    Calculates mass action constant for dissociation of carbonic acid using equation from Table 2.2 in Dreybrodt (1988), originally reported in Plummer and Busenberg (1982).

    Parameters
    ----------
    T_K : float
       temperature Kelvin

    Returns
    -------
    K_2 : float
       mass action constant for dissociation of carbonic acid
    """
    K_1 = 10.**(-356.3094 - 0.06091964*T_K + 21834.37/T_K + 126.8339*np.log10(T_K) - 1684915.0/T_K/T_K)
    return K_1

def calc_K_0(T_K):
    """
    Calculates mass action constant for conversion of CO2 to carbonic acid.

    Calculates mass action constant for conversion of CO2 to carbonic acid using equation from Table 2.2 in Dreybrodt (1988), originally reported in Plummer and Busenberg (1982).

    Parameters
    ----------
    T_K : float
       temperature Kelvin

    Returns
    -------
    K_0 : float
       mass action constant for conversion of CO2 to carbonic acid
    """
    K_1 = 10.**(-356.3094 - 0.06091964*T_K + 21834.37/T_K + 126.8339*np.log10(T_K) - 1684915.0/T_K/T_K)
    K_0 = 1.7*0.0001/K_1
    return K_0

def calc_K_H(T_K):
    """
    Calculates Henry's law constant for CO2.

    Calculates Henry's law constant for CO2 using equation from Table 2.2 in Dreybrodt (1988), originally reported in Plummer and Busenberg (1982).

    Parameters
    ----------
    T_K : float
       temperature Kelvin

    Returns
    -------
    K_H : float
       Henry's law constant for CO2.
    """
    K_H = 10.**(108.3865 + 0.01985076*T_K - 6919.53/T_K - 40.45154*np.log10(T_K) + 669365./T_K/T_K)
    return K_H

def calc_K_W(T_K):
    """
    Calculates mass action constant for dissociation water.

    Calculates mass action constant for dissociation of water using equation from Table 2.2 in Dreybrodt (1988), originally Harned and Hamer (1933).

    Parameters
    ----------
    T_K : float
       temperature Kelvin

    Returns
    -------
    K_W : float
       mass action constant for dissociation of water
    """

    K_W = 10.**(22.801 - 4787.3/T_K - 0.010365*T_K - 7.1321*np.log10(T_K))
    return K_W

def calc_k1(T_K):
    """
    Calculates k1+ kinetic constant from Table 1 of Kaufmann and Dreybrodt (2007). 

    Parameters
    ----------
    T_K : float
       temperature Kelvin

    Returns
    -------
    k1 : float
       kinetic constant k1+

    Notes
    -----
    Related to the rate of CO2 conversion.
    """
    k1 = 10.**(329.850 - 110.54*np.log10(T_K) - 17265.4/T_K)
#    k1 = (10.**-3) * np.exp(34.69 - 9252./T_K)
    return k1

def calc_k2(T_K):
    """
    Calculates k2+ kinetic constant from Table 1 of Kaufmann and Dreybrodt (2007).

    Parameters
    ----------
    T_K : float
       temperature Kelvin

    Returns
    -------
    k2 : float
       kinetic constant k2+

    Notes
    -----
    Related to the rate of CO2 conversion.
    """
#    k2 = 10.**(14.072 - 3025./T_K)
    k2 = 10.**(13.635 - 2895./T_K)
    return k2

def calc_k_neg1(T_K):
    """
    Calculates k1- kinetic constant from Table 1 of Kaufmann and Dreybrodt (2007).

    Parameters
    ----------
    T_K : float
       temperature Kelvin

    Returns
    -------
    k_neg1 : float
       kinetic constant k1-

    Notes
    -----
    Related to the rate of CO2 conversion.
    """
    k_neg1 = 10.**(13.558 - 3617.1/T_K)
    return k_neg1

def calc_k_neg2(T_K):
    """
    Calculates k2- kinetic constant from Table 1 in Kaufmann and Dreybrodt (2007).

    Parameters
    ----------
    T_K : float
       temperature Kelvin

    Returns
    -------
    k_neg2 : float
       kinetic constant k2-

    Notes
    -----
    Related to the rate of CO2 conversion.
    """
    k_neg2 = 10.**(14.09 - 5308./T_K)
    return k_neg2


