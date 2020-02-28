import numpy as np

#properties of common ions
properties = {
    'H':{'weight':1.00794, 'charge':1., 'radius':9.},
    'Ca':{'weight':40.08, 'charge':2., 'radius':6. },
    'Cl':{'weight':35.46, 'charge':-1., 'radius':3.},
    'CO3':{'weight':60.01, 'charge':-2., 'radius':4.5},
    'HCO3':{'weight':61.02, 'charge':-1., 'radius':4.},
    'K':{'weight':39.10, 'charge':1.},
    'Mg':{'weight':24.32, 'charge':2., 'radius':8.},
    'Na':{'weight':22.99, 'charge':1., 'radius':4.},
    'NO3':{'weight':62.01, 'charge':-1.},
    'PO4':{'weight':94.98, 'charge':-3.},
    'SO4':{'weight':96.06, 'charge':-2., 'radius':4.},
    'H2CO3':{'weight':62.03, 'charge':0.},
    'OH':{'weight':17.01, 'charge':-1., 'radius':3.},
    'CO2':{'weight':44.01, 'charge':0},
    'H2CO3s':{'weight':62.03, 'charge':0},
    'CaCO3':{'weight':100.09, 'charge':0}
}

def getProperties():
    """

Returns the dictionary containing ion properties.

"""
    return properties

def condTo25(cond, temp):
    """
    Converts given value of specific conductance to the value at 25 C.

    Uses the relation from Sorensen and Glass (1987), as found in
    Equation 4 of Hayashi (2004).

    Parameters
    ----------
    cond : float or array
       measured value of specific conductance
    temp : float or array
       temperature of water at which specific conductance was measured.
   
    Returns
    -------
    cond25: float or array
       specific conductance at 25 C.

    """
    a = 0.0187 #from Hayashi, 2004
    cond25 = cond/(1+a*(temp - 25.0))
    return cond25

def DebyeHuckel(I, z, r, T = 25.):
    """

    Calculates activity coefficient, gamma, using the Debye-Huckel equation.

    Uses Equation 2.8 in Dreybrodt (1998) (Processes in Karst Systems).

    Parameters
    ----------
    I : float or array
       ionic strength (mol/L)
    z : float or array
       ion charge
    r : float or array
       ionic radius (angstroms)
    T : float or array, optional
       temperature in C, defaults to room temp

    Returns
    -------
    gamma : float or array
       activity coefficient

    """
    A = 0.4883 + 8.074 * T * 10.**(-4)
    B = 0.3241 + 1.6 * T * 10.**(-4)
    log_gamma = -A*(z**2.) * np.sqrt(I) / (1 + B*r*np.sqrt(I))
    gamma = 10**(log_gamma)
    return gamma

#Calculate activity coefficient for a neutral species
# from eq 2.9, Dreybrodt 1988
def neutralGamma(I):
    """
    Calculates the activity coefficient for neutral species.

    Uses Equation 2.9 from Dreybrodt (1988) (Processes in Karst Systems)..

    Parameters
    ----------
    I : float or array
       ionic strength of the solution
    
    Returns
    -------
    gamma : float or array
       ion activity coefficient
    """
    gamma = 10.**(0.1*I)
    return gamma

def gamma(ion, I, T_C=25.):
    """

    Calculates activity coefficient, gamma, for a given ion name.

    Parameters
    ----------
    ion : string
       String containing name of ion.
    I : float
       ionic strength of solution
    T_C : float, optional
       temperature Celsius, default is 25 C

    Returns
    -------
    gamma : float
       ion activity coefficient

    Notes
    -----
    Uses function DebyeHuckel or neutralGamma as appropriate.

    """
    if ion in properties:
        if 'radius' in properties[ion]:
            radius = properties[ion]['radius']
        else:
            if ('charge' in properties[ion]): 
                if (properties[ion]['charge']==0):
                    activity_coef = neutralGamma(I)
                    return activity_coef
            print("Radius of " + ion + " not available.")
            return None
        if 'charge' in properties[ion]:
            charge = properties[ion]['charge']
        else:
            print("Charge of " + ion + " not available.")
            return None
        activity_coef = DebyeHuckel(I,charge,radius,T=T_C)
        return activity_coef
    else:
        print("Ion not in properties list: " + ion)
        return None
    
def HardnessFromCond(cond, T_C=None):
    """

    Function to estimate total hardness as mg/L CaCO3 from solution
    specific conductance (microS/cm) using the empirical equation from
    Krawczyk and Ford (2006), Earth Surface Processes and Landforms.

    Parameters
    ----------
    cond : float or array
       specific conductance in microS/cm (can be float or array of floats)
    T_C : float or array
       temperature in Celsius (can be float or array of floats)

    Returns
    -------
    hardness : float or array
       hardness in mg/L CaCO3 (if cond is an array, an array is returned)

    """
    if np.size(T_C) == 1: #This extra logic allows T_C to be an array
        if T_C != None:
            cond = condTo25(cond, T_C)
    else:
        cond = condTo25(cond, T_C) #If T_C has more than one element
    hardness = (cond - 31.5)/1.86
    return hardness

def CaFromCond(cond, T_C=None, mol_L=False):
    """
    Function to estimate Ca mg/L from solution specific conductance
    (microS/cm) using the empirical equation from Krawczyk and Ford,
    (2006), Earth Surface Processes and Landforms.

    Parameters
    ----------
    cond : float or array 
       specific conductance (microS/cm)
    T_C : float or array 
       temperature in degrees Celsius
    mol_L : bool, optional
      if true then convert from mg/L to mol/L (default false)

    Returns
    -------
    Ca : float or array
       concentration of Ca in mg/L or mol/L
    """
    hardness = HardnessFromCond(cond, T_C=T_C)
    mw_Ca = properties['Ca']['weight']
    mw_CaCO3 = mw_Ca + properties['CO3']['weight']
    Ca = hardness*(mw_Ca/mw_CaCO3)
    if mol_L:
        #convert to mol/L
        Ca = Ca / (properties['Ca']['weight']*1000.)
    return Ca

def CtoK(T_C):
    """Converts Celsius to Kelvin."""
    return T_C + 273.15

def KtoC(T_K):
    """Converts Kelvin to Celsius."""
    return T_K - 273.15

def approxI(metals_conc):
    """
    Calculates approximate ionic strength from the total concentration
    of metals.

    Uses equation 2.19 or 2.20 in Dreybrodt (1988).  Approximation is
    good in most natural karst waters.

    Parameters
    ----------
    metals_conc : float or array
       concentration of metals (mol/L)

    Returns
    -------
    I : float or array
       ionic strength
    """
    I = 3.0 * metals_conc
    return I

def molL_to_mgL(molL, ion):
    """
    Converts concentrations from mol/L to mg/L.
    
    Parameters
    ----------
    molL : float or array
       concentration in mol/L
    ion : str
       name of ion

    Returns
    -------
    mgL : float or array
       concentration in mg/L
    """
    if ion in properties:
        mw = properties[ion]['weight']
        mgL = (1000. * mw)*molL
        return mgL
    else:
        print("Ion +", ion, " not in properties.")
        return -1

def mgL_to_molL(mgL, ion):
    """
    Converts concentrations from mg/L to mol/L.
    
    Parameters
    ----------
    mgL : float or array
       concentration in mg/L
    ion : str
       name of ion

    Returns
    -------
    molL : float or array
       concentration in mol/L
    """
    if ion in properties:
        mw = properties[ion]['weight']
        molL = mgL/(1000.*mw)
        return molL
    else:
        print("Ion +", ion, " not in properties.")
        return -1

def mmolL_to_meqL(mmolL, ion):
    """
    Calculates milliequivalents per liter from mmol concentration.
    
    Parameters
    ----------
    mmolL : float or array
       concentration in mmol/L
    ion : str
       name of ion

    Returns
    -------
    meqL : float or array
       milliequivalents per liter
    """
    if ion in properties:
        charge = properties[ion]['charge']
        meqL = mmolL*charge
        return meqL
    else:
        print("Ion +", ion, " not in properties.")
        return -1

def molL_to_meqL(molL, ion):
    """
    Calculates milliequivalents per liter from molar concentration.
    
    Parameters
    ----------
    molL : float or array
       concentration in mol/L
    ion : str
       name of ion

    Returns
    -------
    meqL : float or array
       milliequivalents per liter
    """
    mmolL = molL*1000.
    meqL = mmolL_to_meqL(mmolL,ion)
    return meqL

##The solution object, which contains solution properties and methods
## for calculations concerning solutions
class solution(object):
    
    """

    A solution class, that is used for calculations related to chemical solutions.

    Parameters
    ----------
    constituents : array like
       a list (or array) of strings with ion names 
    concentrations : array like
       a list (or array) of ion concentrations 
    units : string or list
       a string or list of strings containing concentration units ('mg/L', 'mol/L',  or 'mmol/L')
    T : float
       temperature in degrees Celsius (default=25)
    T_units : str
       temperature units 'C' or 'K' (default='C')
    cond : float
       specific conductance (default=None)
    pH : float 
       pH of solution (default=None)

    Attributes
    ----------
    ions : dict
       A dictionary containing ion activities and concentrations.  The values for a given ion are accessed using a the name of the ion as a key.  Then, the additional keys 'activity', 'conc_mg', 'conc_mol', 'conc_mmol', and 'meq' access the ion activity (in units of mol/L), concentration in mg/L, concentration in mol/L, concentration in mmol/L and meq/L, respectively.
    T_C : float
       temperature of solution in degrees C
    T_K : float
       temperature of solution in degrees K

    """

    def __init__(self, constituents, concentrations, units, T=25., T_units='C', cond=None, pH=None):  
#        print 'units='+str(units)
        if type(units) == str:
            units = [units] * len(constituents)
        if not(len(constituents) == len(units)):            
            if type(units) == list:
                if len(units) == 1:
                    units = [units[0]] * len(constituents)
                else:
                    print("Units must either contain the same number of items as constituents or only contain one item.")
            else:
                print("Units must either be a list or a string.")
        if len(constituents) == len(concentrations):
            self.ions = {}
            for i, ion in enumerate(constituents):
                ion_dict = {}
                if units[i] == 'mg/L': 
                    ion_dict['conc_mg'] = concentrations[i]
                    ion_dict['conc_mol'] = mgL_to_molL(concentrations[i], ion)
                    ion_dict['conc_mmol'] = ion_dict['conc_mol']*1000.
                    ion_dict['meq'] = mmolL_to_meqL(ion_dict['conc_mmol'], ion)
                    self.ions[ion] = ion_dict
                elif units[i] == 'mol/L':
                    ion_dict['conc_mol'] = concentrations[i]
                    ion_dict['conc_mmol'] = 1000.*concentrations[i]                       
                    ion_dict['conc_mg'] = molL_to_mgL(concentrations[i], ion)
                    ion_dict['meq'] = mmolL_to_meqL(ion_dict['conc_mmol'], ion)
                    self.ions[ion] = ion_dict
                elif units[i] == 'mmol/L':
                    ion_dict['conc_mmol'] = concentrations[i]
                    ion_dict['conc_mol'] = concentrations[i]/1000.
                    ion_dict['conc_mg'] = molL_to_mgL(ion_dict['conc_mol'])
                    ion_dict['meq'] = mmolL_to_meqL(ion_dict['conc_mmol'], ion)
                    self.ions[ion] = ion_dict
                else:
                    print("Units for " + ion + " are incorrect, must either be mg/L, mmol/L, or mol/L")
            if T_units == 'C':
                self.T_C = np.float(T)
                self.T_K = CtoK(self.T_C)
            elif T_units == 'K':
                self.T_K = np.float(T)
                self.T_C = KtoC(self.T_K)
            else:
                print("T_units must equal either C or K.")
            if cond != None:
                self.cond = np.float(cond)
            if pH != None:
                self.pH = np.float(pH)
            #loop through ions and calculate activities
            for key, ion_dict in self.ions.items():
                if ('conc_mol' in ion_dict):
                    gamma = self.activity_coef(key) #estimates ionic strength
                    if gamma != None:
                        activity = gamma*ion_dict['conc_mol']
                        self.ions[key]['activity'] = activity
        else:
            print("Length of constituents and concentrations lists must be equal.")

    def I(self):
        """
        Function to calculate ionic strength of the solution object from full equation using known ion properties.

        Returns
        -------
        I - float
           ionic strength
        """
        I=0.
        for ion in self.ions:
            if ion in properties:
                if ('charge' in properties[ion]):
                    I += 0.5*(properties[ion]['charge']**2.)*\
                        self.ions[ion]['conc_mol']
        return I

    def approxI(self):
        """
        Function to calculate approximate ionic strength of solution metal concentration.

        Returns
        -------
        I : float
           ionic strength
        """
        metals_conc = 0.0
        metals = ['Ca', 'Mg']  #could add more metals here later
        for metal in metals:
            if metal in self.ions:
                if 'conc_mol' in self.ions[metal]:
                    metals_conc += self.ions[metal]['conc_mol']
                else:
                    print("No known " + metal + " concentration in mol/L.")
        I = 3.0 * metals_conc
        return I

    def activity(self, ion):
        """
        Retrieve the activity of a given ion in the solution in mol/L.

        Parameters
        ----------
        ion : str
           string containing ion name for which to calculate activity coefficient
        Returns
        -------
        activity : float
           Activity of the specified ion in mol/L. 
        """           
        if (ion in self.ions):
            activity = self.ions[ion]['activity']
            return activity

    def mol(self, ion):
        """
        Retrieve the molar concentration of a given ion in the solution in mol/L.

        Parameters
        ----------
        ion : str
           string containing ion name for which to calculate activity coefficient
        Returns
        -------
        molL : float
           Concentration of the specified ion in mol/L. 
        """           
        if (ion in self.ions):
            conc = self.ions[ion]['conc_mol']
            return conc

    def mmol(self, ion):
        """
        Retrieve the molar concentration of a given ion in the solution in mol/L.

        Parameters
        ----------
        ion : str
           string containing ion name for which to calculate activity coefficient
        Returns
        -------
        mmolL : float
           Concentration of the specified ion in mol/L. 
        """           
        if (ion in self.ions):
            conc = self.ions[ion]['conc_mmol']
            return conc

    def mg(self, ion):
        """
        Retrieve the concentration of a given ion in the solution in mg/L.

        Parameters
        ----------
        ion : str
           string containing ion name for which to calculate activity coefficient
        Returns
        -------
        mgL : float
           Concentration of the specified ion in mg/L. 
        """           
        if (ion in self.ions):
            conc = self.ions[ion]['conc_mg']
            return conc

    def meq(self, ion):
        """
        Retrieve meq for a given ion in the solution.

        Parameters
        ----------
        ion : str
           string containing ion name for which to calculate activity coefficient
        Returns
        -------
        meqL : float
           Concentration of the specified ion in meq/L. 
        """           
        if (ion in self.ions):
            meq = self.ions[ion]['meq']
            return meq
     
    def activity_coef(self, ion, I=None):
        """
        Calculate activity coefficient, gamma, for a given ion in the solution.

        Parameters
        ----------
        ion : str
           string containing ion name for which to calculate activity coefficient
        I : float 
           Ionic strength of solution.  If None, it will be approximated from metal concentration. (default=None)
        Returns
        -------
        gamma : float
           Activity coefficient for the ion specified.
        """
        if I == None:
            #if ionic strength was not input, approximate it
            #for calcite solutions
            I = self.approxI()
        if ion in properties:
            if 'radius' in properties[ion]:
                radius = properties[ion]['radius']
            else:
                if ('charge' in properties[ion]): 
                    if (properties[ion]['charge']==0):
                        activity_coef = neutralGamma(I)
                        return activity_coef
                print("Radius of " + ion + " not available.")
                return None
            if 'charge' in properties[ion]:
                charge = properties[ion]['charge']
            else:
                print("Charge of " + ion + " not available.")
                return None
            activity_coef = DebyeHuckel(I,charge,radius,T=self.T_C)
            return activity_coef
        else:
            print("Ion not in properties list:" + ion)
            return None

    def charge_balance(self, units='%'):
        """
        Calculate the charge balance error for this solution.

        Parameters
        ----------
        units : str (optional) 
           If set to '%' (default) then the function will return the percent charge balance error.  If set to 'meq', then the function will return the difference in cation and anion balance in meq/L.
        Returns
        -------
        balance : float
           The charge balance error for the solution.
        """
        numerator = 0.
        denominator = 0.
        for ion in self.ions:
            numerator = numerator + self.ions[ion]['meq']
            denominator = denominator + abs(self.ions[ion]['meq'])
        if units=='meq':
            return numerator
        elif units=='%':
            return numerator/denominator*100.
        else:
            print("Invalid value for units: ", units)
            return
        
