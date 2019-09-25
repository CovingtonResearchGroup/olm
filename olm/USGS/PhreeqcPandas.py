#! /usr/bin/env python

"""
Collection of functions to interface between WaterChem and PHREEQC.
"""
from __future__ import print_function
from numpy import isnan,size
from scipy.optimize import brentq
from pandas import DataFrame
import pickle as pickle
import os, subprocess, xlrd, sys
from olm.USGS.loadWaterQualityData import loadSiteListData
from olm.general import molL_to_mgL

default_phreeqc_to_WQX_translation = {
    'Temperature, water':'temp',
    'pH':'pH',
    'Alkalinity, total':'Alkalinity',
    'Calcium':'Ca',
    'Magnesium':'Mg',
    'Potassium':'K',
    'Sodium':'Na',
    'Sulfate':'S',
    'Nitrate':'N',
    'Chloride':'Cl'
    }

def readPhreeqcOutput(phreeqcOutputFile):
    """
    Reads the output file from a PHREEQC simulation run.

    Parameters
    ----------
    phreeqcOutputFile : string
       File containing PHREEQC output that should be read.

    Returns
    -------
    simulationDict : dict
       A dictionary containing outputs from the PHREEQC simulation.

    Notes
    -----
    While this function can be used on its own. It is primarily designed to be used by the processPanel function, which provides a more convenient interface.
    """
    with open(phreeqcOutputFile, 'r', encoding='latin_1') as phreeqc_output:
        simulationDict = {}
        #loop through all lines of file
        block = 'beginning'
        endOfRun = False
        for line in phreeqc_output:
            if not(endOfRun):
                if ('ERROR' in line):
                    print("PHREEQC has encountered an error while running this sample.")
                    return None
                if (block == 'beginning'):
#                    print 'beginning'
                    if ('Reading input data for simulation' in line):
                        block = 'simulation'
                elif (block == 'simulation'):
                     if ('Solution composition' in line):
                        block = 'composition'
                     if ('End of run' in line):
                        endOfRun = True
                elif (block == 'composition'):
                    if not('Elements' in line):
                        linestrip = line.strip()
                        linesplit = linestrip.split()
                    if ('Description of solution' in line):
                        block = 'description'
                    elif (len(linesplit) == 3):
                        simulationDict[linesplit[0]] = linesplit[1]
                elif (block == 'description'):
                    if ('Electrical balance' in line):
                        linestrip = line.strip()
                        linesplit = linestrip.split()
                        balance = linesplit[4]
                        simulationDict['Electrical balance'] = balance
                    if ('Percent error' in line):
                        linestrip = line.strip()
                        linesplit = linestrip.split()
                        error = linesplit[4]
                        simulationDict['Percent error'] = error
                    if ('Distribution of species' in line):
                        block = 'distribution'
                elif (block == 'distribution'):
                    if not('Species' in line):
                        linestrip = line.strip()
                        linesplit = linestrip.split()
                        if ('Saturation indices' in line):
                            #this line should not be included in molality, hence the elif
                            block = 'SI'
                        elif (len(linesplit) > 2):
#                            if len(linesplit) == 2:
                                #For now we won't get category totals
                                #simulationDict[linesplit[0]] = linesplit[1]
                            #line to retreive individual species
                            simulationDict[linesplit[0]+'_Molality'] = linesplit[1]
                            simulationDict[linesplit[0]+'_Activity'] = linesplit[2]
#Old code from when we weren't recording individual species but only category totals
#                            else: #grab CO2
#                                if (linesplit[0] == 'CO2'):
#                                    simulationDict[linesplit[0]] = linesplit[1]
                elif (block == 'SI'):
                    if not('Phase' in line):
                        linestrip = line.strip()
                        linesplit = linestrip.split()
                        if (len(linesplit) == 5):
                            simulationDict['SI_'+linesplit[0]] = linesplit[1]
                        if ('End of simulation' in line):
                            #add simulation to list
#                            simData = {'error':simulationError, 'molality':simulationMolality, 'species':simulationSpecies,'SI':simulationSI}
#                            simulationList.append(simData)
                            block = 'simulation'
        return simulationDict


def writePhreeqcInput(sample_row, phreeqc_file_name, phreeqcDict=default_phreeqc_to_WQX_translation, datetext='', charge=None):
    """
    Writes a PHREEQC input file using the row from a site panel.

    Parameters
    ----------
    sample_row : pandas data frame row
       The row from the pandas data frame for the site and date to be processed.

    phreeqc_file_name : str
       Name of PHREEQC input file to write.

    phreeqcDict : dict
       a dictionary with WQX characteristics as keys and phreeqc chemical names as entries. By default, processPanel will use the built in translation dict, default_phreeqc_to_WQX_translation.

    datetext : str
       String containing the text that describes the date as it should be written into the PHREEQC input file.

    charge : str
        String containing name of element (or pH) that should be adjusted to obtain charge balance. This is done internally by PHREEQC. (Default=None)
    Returns
    -------
    status : int
       OK if equal to 1. Error writing to file if equal to -1.
    """
    try:
        phreeqc_input_file = open(phreeqc_file_name, 'w')
        print('SOLUTION ' + ' on ' + datetext, file=phreeqc_input_file)
        print('units mg/l', file=phreeqc_input_file)
        #loop through all characteristics that should be included in PHREEQC analysis
        for characteristic in phreeqcDict.keys():
            #Check to see if this characteristic is in our panel
            if characteristic in list(sample_row.keys()):
                #Check for cases with duplicate dates and no time information
                #In this case, get only first sample on that date
                if size(sample_row.shape)>1 and min(sample_row.shape)>1:
                    sample_row = sample_row.ix[0]
                if not isnan(sample_row[characteristic] ):
                    if phreeqcDict[characteristic]==charge:
                        print(phreeqcDict[characteristic] +' ' + str(sample_row[characteristic])+ ' '+'charge', file=phreeqc_input_file)
                    else:
                        print(phreeqcDict[characteristic] +' ' + str(sample_row[characteristic]), file=phreeqc_input_file)
        print("END", file=phreeqc_input_file)
        phreeqc_input_file.close()
        return 1
    except IOError:
        phreeqc_input_file.close()
        print ("Problem opening sample PHREEQC input file.")
        return -1

def phreeqcRunSetPCO2(logPCO2, phreeqcInputFile, PHREEQC_PATH, DATABASE_FILE, newInputFile=None):
    """
    Reads a PHREEQC input file, modifies it to a set PCO2 and runs the modified file.

    Parameters
    ----------
    phreecInputFile : str
       The name of the input file to modify.

    logPCO2 : float
        Log10 of the PCO2 value to set the sample to.

    PHREEQC_PATH : string
       The path to the PHREEQC executable.

    DATABASE_FILE : string
       The database file to be used by PHREEQC.

    newInputFile : string
        The name of the set-PCO2 file to create and run. If not specified, the final 5 characters will be stripped off of phreeqcInputFile and 'set-PCO2.phrq will be added.'
    Returns
    -------
    simulationDict : Dict
        A dictionary that contains the output of readPhreeqcOutput() run on the set PCO2 simulation.
    """
    with open(phreeqcInputFile) as f:
        input_file_text = f.readlines()
    wrote_CO2_line = False
    PCO2_str = '\tC\t1\tCO2(g)\t' + str(logPCO2)+'\n'
    for linenum, input_line in enumerate(input_file_text):
        #If there is a pH input line we will overwrite it
        if 'pH' in input_line:
            if not wrote_CO2_line:
                input_file_text[linenum] = PCO2_str
                wrote_CO2_line = True
            else:
                print("There is more than one line with pH. This code only works for PHREEQC inputs with a single solution. Check input file!")
    #If there was no pH line, we still need to add the PCO2 line
    if not wrote_CO2_line:
        for linenum, input_line in enumerate(input_file_text):
            if 'SOLUTION' in input_line:
                if not wrote_CO2_line:
                    input_file_text.insert(linenum+1, PCO2_str)
                    wrote_CO2_line = True
                else:
                    print("There is more than one line with SOLUTION. This code only works for PHREEQC inputs with a single solution. Check input file!")
    #Write out new input file text with set PCO2
    if newInputFile == None:
        newInputFile = phreeqcInputFile[:-5]+'-set-PCO2.phrq'
    with open(newInputFile, 'w') as f:
        for line in input_file_text:
            f.write(line)
    #run phreeqc on the input file just created
    blank = ''
    split_file_name = newInputFile.split('.')
    split_list = split_file_name[:-1]#Trim off file ending
    phreeqcOutputFile = blank.join(split_list) + '.out'
    retcode = subprocess.call([os.path.join(PHREEQC_PATH,'phreeqc'), newInputFile, phreeqcOutputFile, DATABASE_FILE])
    print("PHREEQC return code = ",retcode)
    if retcode == 0:
        simulationDict = readPhreeqcOutput(phreeqcOutputFile)
    else:
        #This catches cases where PHREEQC fails (e.g. doesn't converge)
        simulationDict = None
    return simulationDict


    #solution_inputs = readPhreeqcOutput(phreeqcOutputFile)
    #solution_inputs['C']= '   1   CO2(g)   ' + str(logPCO2)
    #writePhreeqcInput(solution_inputs)

def calciteSaturationAtFixedPCO2(logPCO2, phreeqcInputFile, PHREEQC_PATH, DATABASE_FILE, newInputFile=None):
    """
    Function used in root finding of saturation PCO2.

    Function is used by findPCO2atCalciteSaturation(). As a stand alone function, it's better to use phreeqcRunSetPCO2().

    Parameters
    ----------
    phreecInputFile : str
       The name of the input file to modify.

    logPCO2 : float
        Log10 of the PCO2 value to set the sample to.

    PHREEQC_PATH : string
       The path to the PHREEQC executable.

    DATABASE_FILE : string
       The database file to be used by PHREEQC.

    newInputFile : string
        The name of the set-PCO2 file to create and run. If not specified, the final 5 characters will be stripped off of phreeqcInputFile and 'set-PCO2.phrq will be added.'
    Returns
    -------
    SI : float
        Saturation index of calcite for given PCO2.
    """
    simulationDict = phreeqcRunSetPCO2(logPCO2, phreeqcInputFile, PHREEQC_PATH, DATABASE_FILE, newInputFile=newInputFile)
    return float(simulationDict['SI_Calcite'])


def findPCO2atCalciteSaturation(phreeqcInputFile, PHREEQC_PATH, DATABASE_FILE, newInputFile=None, units='atm'):
    """
    Finds the PCO2 where calcite would be saturated by modifing PCO2 of a PHREEQC input file.

    Parameters
    ----------
    phreecInputFile : str
       The name of the input file to modify.

    PHREEQC_PATH : string
       The path to the PHREEQC executable.

    DATABASE_FILE : string
       The database file to be used by PHREEQC.

    newInputFile : string
        The name of the set-PCO2 file to create and run. If not specified, the final 5 characters will be stripped off of phreeqcInputFile and 'set-PCO2.phrq will be added.'

    units: string
        The units in which to return the partial pressure of CO2. Currently, 'ppm' and 'atm' are allowed (Default=atm).
    Returns
    -------
    PCO2 : float
        The PCO2 at which calcite would be in equilibrium.

    """
    args = (phreeqcInputFile, PHREEQC_PATH, DATABASE_FILE, newInputFile)
    eqLogPCO2 = brentq(calciteSaturationAtFixedPCO2, -6., 0., args=args, xtol=0.01)
    if units == 'ppm':
        PCO2 = 1e6 * 10.**eqLogPCO2
    elif units == 'atm':
        PCO2 = 10.**eqLogPCO2
    else:
        print("Specified units were not valid, use 'ppm' or 'atm'.")
        PCO2 = None
    return PCO2

def processPanel(site_panel, site_dir, PHREEQC_PATH, DATABASE_FILE, phreeqcDict=None, force_balance=''):
    """
    Takes a WQXtoPandas site panel and runs all samples through PHREEQC, returning a dataframe of the outputs that is indexed by date. Will run automatically within WQXtoPandas if specified in the excel start file.  However, the function can also be called later by reading in a pickled site panel.

    Parameters
    ----------
    site_panel : Pandas panel
       The Pandas panel object produced for each site by WQXtoPandas that is to be processed through PHREEQC.

    site_dir : string
       The directory containing the site data and where the results should be written out.

    PHREEQC_PATH : string
       The path to the PHREEQC executable.

    DATABASE_FILE : string
       The database file to be used by PHREEQC.

    phreeqcDict : dict
       a dictionary with WQX characteristics as keys and phreeqc chemical names as entries. By default, processPanel will use the built in translation dict, default_phreeqc_to_WQX_translation.

    force_balance : str
       PHREEQC should force charge balance on the ion indicated in this string. Use the string that represents the ion internally in PHREEQC. It is also possible to force balance on pH using force_balance='pH', or on alkalinity using force_balance='Alk'
    """
    if phreeqcDict == None:
        phreeqcDict = default_phreeqc_to_WQX_translation
    output_dates = []
    phreeqc_data_list = []
    df = site_panel['data']
    dates = df.index
    num_converged = 0
    total_num = dates.size
    for date in dates:
        sample_row = df.ix[date]
        datetext = date.date().isoformat()
        if not(os.path.exists(site_dir)):
                    try:
                        os.makedirs(site_dir)
                    except os.error:
                        print(("Problem creating location directory: " + sampledir))
                        return -1
        if force_balance=='':
            phreeqc_file_name = os.path.join(site_dir, datetext+'.phrq')
            writePhreeqcInput(sample_row, phreeqc_file_name, phreeqcDict=phreeqcDict, datetext=datetext)
        else:
            phreeqc_file_name = os.path.join(site_dir, datetext+'-'+force_balance+'.phrq')
            if force_balance=='Alk':
                writePhreeqcInput(sample_row, phreeqc_file_name, phreeqcDict=phreeqcDict, datetext=datetext)
            else:
                writePhreeqcInput(sample_row, phreeqc_file_name, phreeqcDict=phreeqcDict, datetext=datetext, charge=force_balance)
        #run phreeqc on the input file just created
        phreeqcOutputFile = phreeqc_file_name[:-4]+'out'
        retcode = subprocess.call([os.path.join(PHREEQC_PATH,'phreeqc'), phreeqc_file_name, phreeqcOutputFile, DATABASE_FILE])
        print("PHREEQC return code = ",retcode)
        if retcode == 0:
            simulationDict = readPhreeqcOutput(phreeqcOutputFile)
        else:
            #This catches cases where PHREEQC fails (e.g. doesn't converge)
            simulationDict = None
        if force_balance=='Alk' and simulationDict!=None:
            #Force charge balance on Alkalinity
            alk_converged=False
            number_of_iterations = 0
            max_iterations = 10
            while not(alk_converged):
                if (simulationDict != None):
                    #Retrieve charge balance error
                    balance_eq = float(simulationDict['Electrical balance'])
                    print("balance_eq=", balance_eq)
                    if 'Alkalinity' in simulationDict:
                        alk = float(simulationDict['Alkalinity'])
                    else:
                        #Note: Some samples seem to have zero measured alkalinity
                        # might think about conditions to test for under which
                        # convergence doesn't happen or we should quit trying.
                        # Perhaps negative alkalinity??
                        alk = 0.0
                    print("alk=", alk)
                    New_alk_molL = alk + balance_eq#attempt to stabilize this using a factor of 0.9. Otherwise it seems to overshoot. this didn't work.
                    if New_alk_molL < 0:#PHREEQC won't take negative alkalinity inputs
                        New_alk_molL = alk*0.5
                    print("new alk=", New_alk_molL)
                    New_alk_mgL = 0.5*molL_to_mgL(New_alk_molL, 'CaCO3')#PHREEQC wants it as 0.5CaCO3
                    print("new alk mg/L=", New_alk_mgL)
                    sample_row['Alkalinity'] = New_alk_mgL
                    writePhreeqcInput(sample_row, phreeqc_file_name, phreeqcDict=phreeqcDict, datetext=datetext)
                    retcode = subprocess.call([os.path.join(PHREEQC_PATH,'phreeqc'), phreeqc_file_name, phreeqcOutputFile, DATABASE_FILE])
                    print("PHREEQC return code = ",retcode)
                    if retcode == 0:
                        simulationDict = readPhreeqcOutput(phreeqcOutputFile)
                        number_of_iterations = number_of_iterations + 1
                        if abs(float(simulationDict['Percent error']))<5.0:
                            alk_converged=True
                            num_converged = num_converged + 1
                    else:
                        #This catches cases where PHREEQC returns an error (e.g. doesn't converge)
                        simulationDict = None
                    if number_of_iterations>max_iterations:
                        simulationDict=None #This will cause this case to be thrown out below
                else:
                    alk_converged=True
                    print("PHREEQC encountered an error while processing this sample.")
                    phreeqcData = {}
                    phreeqcError = True
                    processThisSample = False
                    reason = "PHREEQC encountered an error while processing this sample. Return code = " + str(retcode)
        if (simulationDict != None):
            #append dates and list of dicts to output lists
            output_dates.append(date) #add this date to index list (row labels)
            phreeqc_data_list.append(simulationDict)
            phreeqcError = False
        else:
            print("PHREEQC encountered an error while processing this sample.")
            phreeqcData = {}
            phreeqcError = True
            processThisSample = False
            reason = "PHREEQC encountered an error while processing this sample. Return code = " + str(retcode)
    #create dataframe from phreeqc outputs
    phreeqc_df = DataFrame(phreeqc_data_list, index=output_dates, dtype='float64')
    if force_balance=='Alk':
        print("Number of converged alkalinity forced charge balance cases = ", num_converged)
        print("Total number of samples = ", total_num)
    return phreeqc_df

def processSites(sitesDir, PHREEQC_PATH, DATABASE_FILE, phreeqcDict=None, regEx='USGS-*', bracket_charge_balance=False, process_regular=True):
    """
    Processes all data from site directories in PHREEQC.

    Parameters
    ----------
    sitesDir : str
       The directory that contains the site directories to be processed.

    PHREEQC_PATH : str
       The path to the phreeqc executable.

    DATABASE_FILE : str
       The path to the phreeqc database file to be used.

    phreeqcDict : dict
       a dictionary with WQX characteristics as keys and phreeqc chemical names as entries. By default, processPanel will use the built in translation dict, default_phreeqc_to_WQX_translation.

    regEx : str (optional)
       A regular expression to be used to locate site directories.  (default = 'USGS-')

    bracket_charge_balance : bool
       If set to True, then charge balance will be bracketed by forcing balance on Ca and Alkalinity alternately. Default = False.
    process_regular : bool
       If set to True, then PHREEQC will be run without any charge balance forcing. This can be set to False if the calculations without forced balance are already done and you only want to do the charge balance runs. Default = True.

    Returns
    -------
    None
    """
    sitesDict = loadSiteListData(regEx=regEx, processedSitesDir = sitesDir)
    if process_regular:
        #Run PHREEQC on sites without forced charge balance
        for site, site_panel in sitesDict.items():
            print(("Processing "+site+" in PHREEQC"))
            sitedf = processPanel(site_panel, os.path.join(sitesDir,site), PHREEQC_PATH=PHREEQC_PATH, DATABASE_FILE=DATABASE_FILE, phreeqcDict=phreeqcDict)
            phreeqc_site_file = os.path.join(sitesDir,site,site+'-PHREEQC.pkl')
            try:
                pickle.dump(sitedf, open(phreeqc_site_file, 'wb'))
                sitedf.to_csv(phreeqc_site_file[:-3]+'csv')
            except IOError:
                print('Problem writing out PHREEQC data file.')
    if bracket_charge_balance:
        #Run PHREEQC using bracketed charge balance for sites
        for site, site_panel in sitesDict.items():
            #Force balance on Calcium
            phreeqc_df_ca = processPanel(site_panel, os.path.join(sitesDir,site), PHREEQC_PATH, DATABASE_FILE, force_balance='Ca')
            phreeqc_site_file_ca = os.path.join(sitesDir,site,site+'-PHREEQC-Ca.pkl')
            try:
                pickle.dump(phreeqc_df_ca, open(phreeqc_site_file_ca, 'wb'))
                phreeqc_df_ca.to_csv(phreeqc_site_file_ca[:-3]+'csv')
            except IOError:
                print('Problem writing out PHREEQC Ca data file.')
            #Force balance on Alkalinity
            phreeqc_df_alk = processPanel(site_panel, os.path.join(sitesDir,site), PHREEQC_PATH, DATABASE_FILE, force_balance='Alk')
            phreeqc_site_file_alk = os.path.join(sitesDir,site,site+'-PHREEQC-Alk.pkl')
            try:
                pickle.dump(phreeqc_df_alk, open(phreeqc_site_file_alk, 'wb'))
                phreeqc_df_alk.to_csv(phreeqc_site_file_alk[:-3]+'csv')
            except IOError:
                print('Problem writing out PHREEQC Alk data file.')


def runPHREEQC(startfilename, process_regular=True):
     """
     Process all samples within a site directory using information from the start file.

     Parameters
     ----------
     startfilename : str
          Name of start file to find settings for running PHREEQC.

     """
     print("Processing site water chemisty data in PHREEQC...")
     startfile = xlrd.open_workbook(startfilename)
     sheet = startfile.sheet_by_index(0)
     #parse start file to determine what should be done
     characteristicsBlockStarted = False
     settingsDict = {}
     for rownum in range(sheet.nrows):
         line = sheet.row_values(rownum)
         if not(characteristicsBlockStarted):
             if not(line[0][0] == '#'): #ignore comments
                 if not(line[0] == 'Characteristic'):
                     settingsDict[line[0]] = line[1]
                 else: #grab the characteristic block column headings
                     characteristicsBlockStarted = True
     sitesDir = os.path.join(
         settingsDict['Path to output directory'],
         settingsDict['Name of output directory'])
     DATABASE_FILE = os.path.join(
         settingsDict['Path to chemical database'],
         settingsDict['Name of chemical database'])
     LOG_FILE = os.path.join(
         settingsDict['Path to output directory'],
         settingsDict['Name of output directory'],
         settingsDict['Log file name'])
     PHREEQC_PATH = settingsDict['Path to PHREEQC']
     bracket_charge_balance = settingsDict['Force balance on Ca and Alk'] == 'Yes'

     processSites(sitesDir, PHREEQC_PATH, DATABASE_FILE, phreeqcDict=None, regEx='USGS-*', bracket_charge_balance=bracket_charge_balance, process_regular=process_regular)


if __name__=="__main__":
    #get site directory  and charge bracketing from command line argument
    startfilename = sys.argv[1]
    process_regular = True
    if len(sys.argv)>2:
        if sys.argv[2]=='False':
            process_regular=False
    runPHREEQC(startfilename, process_regular=process_regular)
