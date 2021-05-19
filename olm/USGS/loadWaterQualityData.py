"""
Functions to load water quality data that has been processed and pickled by WQXtoPHREEQC
"""
import os
import pickle as pickle
from pandas.io.pickle import read_pickle
from .siteListExtraction import *
from glob import glob

DEFAULT_DIR = './Processed-Sites'

def loadSiteListData(siteListText = None,
                     siteFile = None,
                     regEx = 'USGS-*',
                     processedSitesDir = DEFAULT_DIR,
                     loadPhreeqc = False,
                     loadMetaData = False
                     ):
    """
    Retrieves site data for multiple sites within a processed sites directory.

    Parameters
    ----------
    siteListText : string (optional)
        a list of sites separated by semi-colons

    siteFile : string (optional)
        a filename of a text file with a list of sites

    regEx : string (optional)
        regular expression used to search for site directories within the processed sites directory (default = 'USGS-')

    processedSitesDir : string (optional)
        directory that contains all of the processed site directories. It is important to change this if the default is not correct. (default='./Processed-Sites')

    loadPhreeqc : boolean
        If set to true, PHREEQC outputs will also be loaded for each site. (default=False)

    loadMetaData : boolean
        If set to true, the site metadata will be loaded for each site. (default=False)

    Returns
    -------
    sitesDict : dict
        A dictionary of site DataFrames keyed by site name.

    or if loadPhreeqc or loadMetaData are set to true
    (sitesDict, sitesPheeqcDict, sitesMetaDataDict) : tuple
       A tuple containing the sitesDict and dicts of the PHREEQC data and/or metadata for each site. Order is as shown.

    """
    siteList = -1
    #If the needed data is provided to find the site list then use it
    if not(siteListText == None):
        #check whether we have a valid site directory
        processedSitesDir = checkSitesDir(processedSitesDir)
        siteList = siteListFromLine(siteListText)
    elif not(siteFile == None):
        #check whether we have a valid site directory
        processedSitesDir = checkSitesDir(processedSitesDir)
        siteList = siteListFromFile(siteFile)
    elif not(regEx == None):
        #check whether we have a valid site directory
        processedSitesDir = checkSitesDir(processedSitesDir)
        siteList = siteListFromRegEx(regEx, processedSitesDir = processedSitesDir)
    else:
        # if not provided then query user for needed data
        processedSitesInput = input("Path of the processed sites directory (Default = ./Processed-Sites): ")
        if (processedSitesInput != ''):
            processedSitesDir = processedSitesInput
        print(processedSitesDir)
        processedSitesDir = checkSitesDir(processedSitesDir)
        modeOK = False
        while not(modeOK):
            mode = input("Do you want to: \n \t 1) enter a semi-colon separated list of sites \n \t 2) provide a text file of sites \n \t \n\t 3) provide an XML list of sites \n \t 4) provide a wildcard expression to obtain sites from directory list\n Enter 1, 2, 3, or 4: ")
            if mode.isdigit():
                if ( (int(mode) > 0) and (int(mode) < 5) ):
                    modeOK = True
                else:
                    print("Invalid input")
            else:
                print("Invalid input")
        if (int(mode) == 1):
            siteListText = input("Enter list of sites separated by semi-colons: ")
            siteList = siteListFromLine(siteListText)
        elif (int(mode) == 2):
            siteFile = input("Enter path to text file containing site list: ")
            siteList = siteListFromFile(siteFile)
        elif (int(mode) == 3):
            siteFile = input("Enter path to XML file containing site list: ")
            siteList = siteListFromFile(siteFile, XML=True)
        elif (int(mode) == 4):
            regEx = input("Enter regular expression: ")
            siteList = siteListFromRegEx(regEx)
    if (siteList != -1):
        #process the sites in the list

        sitesDict = {}
        sitesPhreeqcDict = {}
        sitesMetaDataDict = {}
        for site in siteList:
            siteFrame = loadSiteData(site, processedSitesDir = processedSitesDir)
            if siteFrame is not None: #If site data does not read in correctly, loadSiteData returns None
                sitesDict[site] = siteFrame
                if loadPhreeqc:
                    sitedf = loadSitePhreeqcData(site, processedSitesDir = processedSitesDir)
                    sitesPhreeqcDict[site] = sitedf
                if loadMetaData:
                    siteMetaData = loadSiteMetaData(site, processedSitesDir = processedSitesDir)
                    sitesMetaDataDict[site] = siteMetaData
        if loadPhreeqc or loadMetaData:
            return_list = [sitesDict]
            if loadPhreeqc:
                return_list.append(sitesPhreeqcDict)
            if loadMetaData:
                return_list.append(sitesMetaDataDict)
            return tuple(return_list)
        else:
            return sitesDict

def loadSiteMetaData(site, processedSitesDir = DEFAULT_DIR):
    #Add USGS tag if needed
#    if not(site.startswith('USGS-')):
#        site = 'USGS-'+site
    try:
        metaDataFile = os.path.join(processedSitesDir, site, site+'-Site-Description.pkl')
        siteMetaData = pickle.load(open(metaDataFile, 'rb'))
    except IOError:
        print(("Problem reading pickle file: " + metaDataFile ))
        return None
    return siteMetaData



def loadSiteData(site, processedSitesDir = DEFAULT_DIR):
    """
    Retrieves site data for an individual site from a directory of processed sites.

    Parameters
    ----------
    site : string
        name of site to retrieve, with or without USGS- tag at beginning.

    processedSitesDir : string (optional)
        directory that contains the processed site directory associated with the desired site. It is important to change this if the default is not correct. (default='./Processed-Sites')

    Returns
    -------
    siteDataFrame : pandas.core.dataframe.DataFrame
        A pandas multiindexed DataFrame object with data and metadata from the requested site.

    """
    #Add USGS tag if needed
#    if not(site.startswith('USGS-')):
#        site = 'USGS-'+site
    try:
        frameFile = os.path.join(processedSitesDir, site, site+'-Dataframe.pkl')
        siteFrame = read_pickle(frameFile)
    except IOError:
        print(("Problem reading pickle file: " + frameFile ))
        return None
    return siteFrame

def loadSitePhreeqcData(site, processedSitesDir = DEFAULT_DIR):
    """
    Retrieves site PHREEQC data for an individual site from a directory of processed sites.

    Parameters
    ----------
    site : string
        name of site to retrieve, with or without USGS- tag at beginning.

    processedSitesDir : string (optional)
        directory that contains the processed site directory associated with the desired site. It is important to change this if the default is not correct. (default='./Processed-Sites')

    Returns
    -------
    sitedf : pandas.core.frame.DataFrame
        A pandas dataframe object with PHREEQC data from the requested site.

    """
    #Add USGS tag if needed
#    if not(site.startswith('USGS-')):
#        site = 'USGS-'+site
    try:
        phreeqcFile = os.path.join(processedSitesDir, site, site+'-PHREEQC.pkl')
        sitedf = read_pickle(phreeqcFile)
    except IOError:
        print(("Problem reading pickle file: " + phreeqcFile ))
        return None
    return sitedf

def siteListFromLine(siteListText):
    siteList = siteListText.split(';')
    siteList = [x.strip() for x in siteList]
    #check for USGS Tag at beginning of site number
#    for i, site in enumerate(siteList):
#        if not(site.startswith('USGS-')):
#            siteList[i] = 'USGS-' + siteList[i]
    return siteList

def siteListFromFile(siteFile,
                     sitesDir = DEFAULT_DIR,
                     XML=False):
    if (siteFile.endswith('.xml') or (XML == True)):
        siteList = extractSitesFromXML(siteFile)
    else:
        siteList = extractSitesFromText(siteFile)
    #check for USGS Tag at beginning of site number
#    for i, site in enumerate(siteList):
#        if not(site.startswith('USGS-')):
#            siteList[i] = 'USGS-' + siteList[i]
 #   siteList = [os.path.join(processedSitesDir, x) for x in siteList]
    return siteList

def siteListFromRegEx(regEx,
                      processedSitesDir = DEFAULT_DIR):
#    print("processedSitesDir="+processedSitesDir)
    listText = os.path.join(processedSitesDir, regEx)
    sitePath = glob(listText)
    siteList = []
    for site in sitePath:
        if os.path.isdir(site):
            head,tail = os.path.split(site)
            siteList.append(tail)
    return siteList

def checkSitesDir(processedSitesDir):
    while not os.path.exists(processedSitesDir):
        print("Invalid path to processed sites directory.")
        processedSitesDir = input("Path of the processed sites directory (Default = ./Processed-Sites): ")
    return processedSitesDir
