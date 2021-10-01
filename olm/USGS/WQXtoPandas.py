#! /usr/bin/env python

"""
This script takes an input file in USGS/EPA WQX xml format and creates a multi-indexed Pandas Dataframe that contains
time series of water quality data and discharge combined with layers that contain meta data for each data value.
You can call the script from the command line using WQXtoPandas [input WQX start file], or import the runWQXtoPandas
function for calling within a Python session.
"""

from __future__ import print_function
import sys,xlrd,os,subprocess,string,requests
from glob import glob
from math import ceil
import pickle as pickle
from lxml import etree
from pandas import DataFrame, to_datetime, Series, concat, ExcelWriter
from olm.USGS.PhreeqcPandas import processPanel

#import functions from olm package
from olm.USGS.siteListExtraction import extractSitesFromXML
from olm.USGS.siteListExtraction import extractSitesFromText
from olm.USGS.DataRetrieval import querySiteList, GetDailyDischarge, GetSiteData
from olm.USGS.dataSlice import extractValues

def WQXtoPandas(xmlLocation, charDict, outputPath='.', fromFile=False, outputDirName='Processed-Sites',
                RUN_PHREEQC=False, PHREEQC_PATH='/home/mcoving/phreeqc-2.18.0/bin/',
                DATABASE_FILE='/home/mcoving/phreeqc-2.18.0/database/phreeqc.dat', LOG_FILE = 'Result.log',
                START_FILE = None, splittag='',bracket_charge_balance=False):
    """
    Processes a WQX xml data file and loads data for each site in the WQX file into Pandas data objects that are
    stored in directories for each site.

    Parameters
    ----------
    xmlLocation : string
       Content depends on mode in which WQXtoPandas is run. When fromFile is set to False (input methods 2 or 3 in excel file) this string contains the html for a query to the USGS NWIS database to obtain an xml file of the desired data.  Alternatively, if fromFile is True (input method 1 in excel file) then this string contains the name of the xml file from which to read the data.

    charDict : dict
       A dictionary containing information about the characteristics to be processed.  Keys are EPA SRS characteristic names. Each entry in the dictionary is a second dictionary that contains keys IsRequired, pcode, fraction, and quality. These entries tell WQXtoPandas whether a given characteristic is required in order to process a sample, and whether a specific pcode, fraction, or quality should be required.  See excel example file for more details.

    outputPath : string
       path to directory that will contain output directory

    fromFile : boolean
       True if data will be read from an xml file already present on computer.  False if xml file should be queried from NWIS. (Default=False)

    outputDirName : string
       Name of output directory where all site data will be written out. (Default='Processed-Sites')

    RUN_PHREEQC : boolean
       Set to true if samples should be processed through PHREEQC. (Default=False)
    PHREEQC_PATH : string
       Path to PHREEQC executable (folder only, not executable file name)

    DATABASE_FILE : string
       Path to database file that PHREEQC should use, including database file name.
    LOG_FILE : string
       Name of log file that WQXtoPandas will create. (Default='Result.log')

    START_FILE : string
       Name of xls start file that was used to run this instance of WQXtoPandas. Name will be written out in log file.

    bracket_charge_balance : bool
       If set to true, WQXtoPandas will alternately force charge balance on calcium and alkalinity, while the latter is not physically meaningful, this provides a useful estimate of uncertainty for cases with high charge balance errors.  This is most useful for water that is very dilute or with high organic content, such that titrated alkalinity values are artificially high.

    Returns
    -------

    Returns 0 if execution successful.  Returns -1 in case of error.

    Notes
    -----

    Designed to be run through convenience function runWQXtoPandas().
    """
    try:
        #Check to see if output directory exists
        absOutputDirPath = os.path.abspath(outputPath)
        sitesdir = os.path.join(absOutputDirPath, outputDirName)
        print("sitesdir",sitesdir)
        if not(os.path.exists(sitesdir)):
            try:
                os.makedirs(sitesdir)
            except os.error:
                print(("Problem creating output directory. Check output path name: "+outputPath))
                return -1
        #create xml tree
        if fromFile:
            #read from file
            wqxtree = etree.ElementTree(file=xmlLocation)
        else:
            #check whether we already have a matching xml file
            xmlSaveFile = LOG_FILE + splittag + '.xml'
            if ( os.path.isfile(xmlSaveFile) ):
                goodAnswer = False
                while not(goodAnswer):
                    answer = input("An xml file ("+ xmlSaveFile+ ") already exists.  \n Use this instead of html query (y or n)?")
                    if (answer.startswith('y')):
                        #read from file
                        wqxtree = etree.ElementTree(file=xmlSaveFile)
                        goodAnswer = True
                        queryXML = False
                    elif (answer.startswith('n')):
                        goodAnswer = True
                        queryXML = True
            else:
                queryXML = True
            #If we don't have a matching xml file, or we want to obtain a new one, then get the new xml
            if (queryXML):
                print("Obtaining xml file from USGS NWIS using html query...")
                #parse from html query
                print("XML query string: ",xmlLocation)
                r = requests.get(xmlLocation)
                if not r.ok:
                    #There is some problem with the xml query
                    print("Response: ", str(r))
                    print("Reason: ", r.reason)
                    print("Warning: ", r.headers['Warning'])
                #write to xml file
                try:
                    #write xml to file
                    xmlFile = open(xmlSaveFile, 'w')
                    print(r.text, file=xmlFile)
                    xmlFile.close()
                    wqxtree = etree.ElementTree(file = xmlSaveFile)
                except IOError:
                    print(("Problem writing to xml file to store html query: " + xmlSaveFile))
                    return -1
        #begin parsing XML tree
        root = wqxtree.getroot()
        #get namespace map
        NSMAP = root.nsmap
        WQX = "{%s}" % NSMAP[None]
        #iterate over all <Activity> tags within file and process each sample
        samples_processed = []
        samples_not_processed = []
        sitesDict = {}
        for activity in wqxtree.getiterator(tag=WQX + "Activity"):
            processThisSample = True
            reason = ''
            description = activity.find(WQX + "ActivityDescription")
            if (description != None):
                datetext = description.findtext(WQX + "ActivityStartDate")
                starttime = description.find(WQX + "ActivityStartTime")
                if (starttime != None):
                    timetext = starttime.findtext(WQX + "Time")
                    timezone = starttime.findtext(WQX + "TimeZoneCode")
                else:
                    timetext = ''
                    timezone = ''
                location = description.findtext(WQX + "MonitoringLocationIdentifier")
                if (location[:5] =='USGS-'):
                    USGS=True
                else:
                    USGS=False
                descriptionDict = {'location':location, 'date':datetext, 'time':timetext, 'timezone':timezone}
            else:
                descriptionDict = None
                processThisSample = False
                reason = 'No description'
            print(('Processing sample from ' + location + ' on ' + datetext))
            #create null sample dict
            sampleDict = {}
            sampleMetaDict = {}
            #iterate though all results for this activity
            for result in activity.getiterator(tag = WQX + 'Result'):
                if (processThisSample):
                    try:
                        resultdesc = result.find(WQX + "ResultDescription")
                        characteristic = resultdesc.findtext(WQX + "CharacteristicName")
                        if (characteristic in charDict):
                            samplefraction = resultdesc.findtext(WQX + "ResultSampleFractionText")
                            pcode = resultdesc.findtext(WQX + "USGSPCode")
                            quality = resultdesc.findtext(WQX + "ResultStatusIdentifier")
                            measure = resultdesc.find(WQX + "ResultMeasure")
                            count = 1.0
                            detection = resultdesc.findtext(WQX + "ResultDetectionConditionText")
                            #print('detection=',detection)
                            if not(measure == None) or not(detection == None):
                                if not(measure == None):
                                    value = measure.findtext(WQX + "ResultMeasureValue")
                                    #print('initial value = ',value)
                                    units = measure.findtext(WQX + "MeasureUnitCode")
                                    #EPA system does not have detection info.
                                    #Check for < in value text.
                                    if '<' in str(value):
                                        value = value[1:]
                                        nondetect = True
                                    else:
                                        nondetect = False
                                elif not(detection == None):
                                    #print("entering nondetect...")
                                    nondetect = True
                                    value = None
                                    labinfo = result.find(WQX + "ResultLabInformation")
                                    if not(labinfo==None):
                                        #print("labinfo present")
                                        quantLimitMeasure = labinfo.find(WQX + "ResultDetectionQuantitationLimit")
                                        if not(quantLimitMeasure==None):
                                            #print("Quant limit present")
                                            nondetectmeasure = quantLimitMeasure.find(WQX + "DetectionQuantitationLimitMeasure")
                                            if not(nondetectmeasure == None):
                                                #print("Quant limit measure present")
                                                value = nondetectmeasure.findtext(WQX + "MeasureValue")
                                                #print('measurevalue=',value)
                                                #print(quantLimitMeasure)
                                    #print("nondetect value=",value)
                                #split pcode into list
                                tempPcodeList = charDict[characteristic]['pcode'].split(';')
    #                            print("tempPcodeList="+str(tempPcodeList))
                                pcodeDict = {}
                                for codePriority, code in enumerate(tempPcodeList):
                                    code = code.strip()
                                    if code != '':
                                        pcodeDict[code] = codePriority
                                #Check whether characteristic meets criteria
                                #for inclusion, otherwise don't add to sampleDict
                                addCharacteristic = True
                                if (charDict[characteristic]['fraction'] != '0'):
                                    #test for correct fraction
                                    if (charDict[characteristic]['fraction'] != samplefraction):
                                        addCharacteristic= False
                                if (addCharacteristic):
                                    if USGS:
                                        if (charDict[characteristic]['pcode'] != '0'):
                                            #test for correct pcode
    #                                        print("pcode = "+pcode)
    #                                        print("pcodeList = "+str(pcodeList))
    #                                        print("pcode in list="+str(pcode in pcodeList))
                                            if not(pcode in pcodeDict):
                                                addCharacteristic = False
                                if (addCharacteristic):
                                    if (charDict[characteristic]['quality'] != '0'):
                                        #test for correct data quality
                                        if (charDict[characteristic]['quality'] != quality):
                                           addCharacteristic = False
                                #end of characteristic criteria check
                                #Process duplicate characteristics
                                if (addCharacteristic):
                                    if (characteristic in sampleDict):
                                        if USGS:
                                            priorPcode =sampleMetaDict[characteristic]['pcode']
                                            #if there are already multiple pcodes get only first one
                                            priorPcode = priorPcode.split(';')[0]
                                            averageValue = False
                                            if (len(pcodeDict) > 1):
                                                thisPcodePriority = pcodeDict[pcode]
                                                priorPcodePriority = \
                                                    pcodeDict[priorPcode]
                                                if (thisPcodePriority >\
                                                        priorPcodePriority):
                                                    #previous characteristic remains
                                                    addCharacteristic = False
                                                elif (thisPcodePriority ==\
                                                      priorPcodePriority):
                                                    averageValue = True
                                            else:
                                                averageValue = True
                                            if averageValue:
                                                #average this value with existing values
                                                count = \
                                                    sampleMetaDict[characteristic]['count']
                                                count += 1.
                                                oldvalue = float(\
                                                    sampleDict[characteristic])
                                                newvalue = (oldvalue * (count - 1.)\
                                                                + float(value))/count
                                                value = str(newvalue)
                                                pcode = priorPcode + '; '+ pcode
                                                priorUnits = \
                                                    sampleMetaDict[characteristic]['units']
                                                units = priorUnits + '; ' + units

                                if (addCharacteristic):
                                    sampleDict[characteristic] = value
                                    sampleMetaDict[characteristic] = {'samplefraction':samplefraction, 'units':units, 'pcode':pcode, 'quality':quality, 'count':count, 'nondetect':nondetect}
                    except etree.XMLSyntaxError as detail:
                        print("File contains invalid XML syntax: ", detail)
                        processThisSample = False
                        reason = "Entry contains invalid XML syntax."
            #end results loop
            #check whether sample has all the required constituents
            if (processThisSample):
                for characteristic in charDict.keys():
                    if (charDict[characteristic]['IsRequired'] != '0'):
                        if not(characteristic in sampleDict):
                            processThisSample = False
                            reason += characteristic+ ' not available. '
            if (processThisSample):
                #check to see whether site directory exists, if not, create it
                sampledir = os.path.join(sitesdir, location)
                if not(os.path.exists(sampledir)):
                    try:
                        os.makedirs(sampledir)
                    except os.error:
                        print(("Problem creating location directory: " + sampledir))
                        processThisSample = False
                        reason = "Problem creating location directory: " + sampledir

            if (processThisSample):
                #Pull daily discharge data from USGS website
                good_discharge_value = False
                num_Q_tries = 0
                if not USGS:
                    #We do not have a USGS site, do not query discharge
                    num_Q_tries = 99
                    dischargeDict = None

                #Try 5 times to retrieve discharge value
                while (not good_discharge_value) and num_Q_tries<=5:
                    dischargeDict = GetDailyDischarge(location, datetext) #currently hard-wired to pcode 00060 (daily discharge, cfs)
                    if dischargeDict != -1:
                        good_discharge_value = True
                    else:
                        num_Q_tries += 1
                        dischargeDict = None
                if (dischargeDict != None):
                    sampleDict['Stream flow, mean. daily'] = dischargeDict['discharge']
                    sampleMetaDict['Stream flow, mean. daily'] = {'units':'cfs', 'pcode':'00060', 'quality':dischargeDict['quality'], 'count':1, 'samplefraction':None, 'nondetect':False}
                    descriptionDict['name'] = dischargeDict['name']
                else:
                    #Possibly allow this sample to be thrown out if no mean daily discharge, and/or similar for instantaneous discharge
                    sampleDict['Stream flow, mean. daily'] = None
                    sampleMetaDict['Stream flow, mean. daily'] = {'units':'cfs', 'pcode':'00060', 'quality':None, 'count':1, 'samplefraction':None, 'nondetect':False}
                # Create data frame row for this sample date
                if descriptionDict['time'] != '':
                    rowdate = to_datetime(datetext+' '+descriptionDict['time'])
                else:
                    rowdate = to_datetime(datetext)
                #Create Multiindex Dataframe to contain sample meta data
                sampleMultiindexRow = concat({
                        'data':DataFrame(sampleDict, index=[rowdate], dtype='float'),
                        'time':DataFrame(descriptionDict['time'],index=[rowdate], columns=list(sampleMetaDict.keys())),
                        'timezone':DataFrame(descriptionDict['timezone'],index=[rowdate], columns=list(sampleMetaDict.keys())),
                        'pcode':DataFrame([extractValues(sampleMetaDict, ['pcode'])['values']], index=[rowdate], columns=list(sampleMetaDict.keys())),
                        'quality':DataFrame([extractValues(sampleMetaDict, ['quality'])['values']], index=[rowdate], columns=list(sampleMetaDict.keys())),
                        'fraction':DataFrame([extractValues(sampleMetaDict, ['samplefraction'])['values']], index=[rowdate], columns=list(sampleMetaDict.keys())),
                        'units':DataFrame([extractValues(sampleMetaDict, ['units'])['values']], index=[rowdate], columns=list(sampleMetaDict.keys())),
                        'count':DataFrame([extractValues(sampleMetaDict, ['count'])['values']], index=[rowdate], columns=list(sampleMetaDict.keys())),
                        'nondetect':DataFrame([extractValues(sampleMetaDict, ['nondetect'])['values']], index=[rowdate], columns=list(sampleMetaDict.keys())),
                        },
                        axis=1)
                #sampleMetaRow = Series(sampleMetaDict, index=[to_datetime(datetext)], dtype='object')
                #Previous solution was reading/writing from pickle files
                #New solution will keep all data in memory until end.
                #This could cause memory problems with large data sets

                #Test whether a df for this location already exists
                if location in sitesDict:
#                    tempDF = sitesDict[location]
#                    sitesDict[location] = tempDF.append(sampleRow)
                    tempMultiindex = sitesDict[location]
                    sitesDict[location] = concat([tempMultiindex, sampleMultiindexRow], axis=0)
                else:
                    sitesDict[location] = sampleMultiindexRow
            #add one to number of samples processed
            if (processThisSample):
                samples_processed.append(location + ' ' + datetext)
            else:
                samples_not_processed.append(location + ' ' + datetext + ' - ' + reason)
        print(('Number of Samples Processed = '+ str(len(samples_processed))))
        print(('Number of Samples Not Processed = ' + str(len(samples_not_processed))))

        #Write out individual site data pickle and csv files in each site directory
        print('Writing out site data files...')
        for location, midf in sitesDict.items():
            print(location)
            pickleFile =  os.path.join(sitesdir, location, location + '-Dataframe.pkl')
            pickle.dump(midf, open(pickleFile, 'wb'))
            #midf.to_excel(pickleFile[:-3]+'xls')
            midx = midf.keys()
            with ExcelWriter(pickleFile[:-3]+'xlsx') as writer:
                for sheet in midx.droplevel(level=1).drop_duplicates().values:
                    midf[sheet].to_excel(writer, sheet_name=sheet)
            #Retrieve and store site description metadata
            siteDescriptionDataDF = GetSiteData(location)
            siteDescriptionDataFileName = os.path.join(sitesdir,location,location+'-Site-Description.pkl')
            pickle.dump(siteDescriptionDataDF, open(siteDescriptionDataFileName, 'wb'))
            siteDescriptionDataDF.to_csv(siteDescriptionDataFileName[:-3]+'csv')
        #Process sites through PHREEQC
        if RUN_PHREEQC:
            print("Processing site water chemisty data in PHREEQC...")
            for location, pnl in sitesDict.items():
                phreeqc_df = processPanel(pnl, os.path.join(sitesdir,location), PHREEQC_PATH, DATABASE_FILE)
                phreeqc_site_file = os.path.join(sitesdir,location,location+'-PHREEQC.pkl')
                try:
                    pickle.dump(phreeqc_df, open(phreeqc_site_file, 'wb'))
                    phreeqc_df.to_csv(phreeqc_site_file[:-3]+'csv')
                except IOError:
                    print('Problem writing out PHREEQC data file.')
            if bracket_charge_balance:
                for location, pnl in sitesDict.items():
                    #Force balance on Calcium
                    phreeqc_df_ca = processPanel(pnl, os.path.join(sitesdir,location), PHREEQC_PATH, DATABASE_FILE, force_balance='Ca')
                    phreeqc_site_file_ca = os.path.join(sitesdir,location,location+'-PHREEQC-Ca.pkl')
                    try:
                        pickle.dump(phreeqc_df_ca, open(phreeqc_site_file_ca, 'wb'))
                        phreeqc_df_ca.to_csv(phreeqc_site_file_ca[:-3]+'csv')
                    except IOError:
                        print('Problem writing out PHREEQC Ca data file.')
                    #Force balance on Alkalinity
                    phreeqc_df_alk = processPanel(pnl, os.path.join(sitesdir,location), PHREEQC_PATH, DATABASE_FILE, force_balance='Alk')
                    phreeqc_site_file_alk = os.path.join(sitesdir,location,location+'-PHREEQC-Alk.pkl')
                    try:
                        pickle.dump(phreeqc_df_alk, open(phreeqc_site_file_alk, 'wb'))
                        phreeqc_df_alk.to_csv(phreeqc_site_file_alk[:-3]+'csv')
                    except IOError:
                        print('Problem writing out PHREEQC Alk data file.')
        #Create log file
        print(('Writing log file: '+LOG_FILE+splittag))
        try:
            log_file = open(LOG_FILE+splittag, 'w')
            print('Start file = ' + START_FILE, file=log_file)
            print('Number of Samples Processed = '+ str(len(samples_processed)), file=log_file)
            print('Number of Samples Not Processed = ' + str(len(samples_not_processed)), file=log_file)
            print("###############", file=log_file)
            print("Characteristics", file=log_file)
            print("###############", file=log_file)
            printColumnNames = True
            for key, flags in charDict.items():
                if (printColumnNames):
                    names = ['characteristic']# + '\t'
                    for column in flags.keys():
                        names.append(str(column))
                    print(str("\t".join(names)), file=log_file)
                    printColumnNames = False
                columns = [key]
                for column in flags.keys():
                    if isinstance(flags[column], str):
                        columns.append(flags[column])
                print(str("\t".join(columns)), file=log_file)
            print("###############", file=log_file)
            print("Samples processed", file=log_file)
            print("###############", file=log_file)
            for line in samples_processed:
                print(line, file=log_file)
            print("###############", file=log_file)
            print("Samples not processed", file=log_file)
            print("###############", file=log_file)
            for line in samples_not_processed:
                print(line, file=log_file)
        except IOError:
            print(("Problem opening log file: " + LOG_FILE))
            return -1
    #exceptions for parsing of xml file
    except IOError:
        print ("Error opening xml file. Does it exist?")
        #Note: can throw this error when discharge values are not read correctly,
        #I should fix this, 6/16/2014
    except etree.XMLSyntaxError as detail:
        print("File contains invalid XML syntax: ", detail)
    except requests.exceptions.RequestException as detail:
        print("Error retrieving data by xml query: ", detail)
    return 0


def runWQXtoPandas(startfilename, autosplitnum=20):
    """
    Runs WQXtoPandas on an excel format input file where parameters can be set for an automatic query of data from
    the USGS NWIS database.

    Parameters
    ----------
    startfilename : string
        A string containing the name of the excel file to be used for input parameters to WQXtoPandas

    autosplitnum : int (optional)
        The number of sites at which a NWIS query is split into multiple queries. (default=20)

    Returns
    -------
    None

    Notes
    -----

    Can be run from within a python shell or script, or as a standalone script from the command line where the start
    file name is provided as the first command line argument (e.g. WQXtoPandas <start file name> <autosplitnum>).
    """
    #PHREEQC input file path
    PHREEQC_INPUT_PATH = './'
    num_samples = 0
    num_processed = 0
    if not(type(autosplitnum)==int):
        print("autosplitnum must be an integer.")
        return -1
    print(('Processing: '+ startfilename))
    try:
        #open start file
        startfile = xlrd.open_workbook(startfilename)
        #open sheet
        sheet = startfile.sheet_by_index(0)
        #parse start file to determine what should be done
        characteristicsBlockStarted = False
        settingsDict = {}
        charDict = {}
        for rownum in range(sheet.nrows):
            line = sheet.row_values(rownum)
            if not(line[0][0] == '#'): #ignore comments
                if not(characteristicsBlockStarted): #read script settings
                    if not(line[0] == 'Characteristic'):
                        settingsDict[line[0]] = line[1]
                    else: #grab the characteristic block column headings
                        column_headings = line[1:]
                        characteristicsBlockStarted = True
                else: #we are in the characteristics block
                    charDict[line[0]] = dict(list(zip(column_headings,line[1:])))
        DATABASE_FILE = os.path.join(
            settingsDict['Path to chemical database'],
            settingsDict['Name of chemical database'])
        LOG_FILE = os.path.join(
            settingsDict['Path to output directory'],
            settingsDict['Name of output directory'],
            settingsDict['Log file name'])
        RUN_PHREEQC = settingsDict['Run PHREEQC?'] == 'Yes'
        bracket_charge_balance = settingsDict['Force balance on Ca and Alk'] == 'Yes'
        if (settingsDict['Input method'] == '1'):
            #We already have an XML file to process that contains water quality data
            #Check whether a wildcard was used and more than one xml file is available
            xml_file_string = os.path.join(
                settingsDict['Path to output directory'],
                settingsDict['Name of output directory'],
                settingsDict['Input file'])
            xml_list = glob(xml_file_string)
            if xml_list==[]:
                print("Empty xml file list. Check path for xml file.")
                print("xml file string =", xml_file_string)
                return -1
            n_xml = len(xml_list)
            if n_xml>1:
                for xml_file in xml_list:
                    WQXtoPandas(
                        xml_file,
                        charDict,
                        outputPath = settingsDict['Path to output directory'],
                        outputDirName = settingsDict['Name of output directory'],
                        fromFile = True,
                        RUN_PHREEQC = RUN_PHREEQC,
                        bracket_charge_balance=bracket_charge_balance,
                        PHREEQC_PATH = settingsDict['Path to PHREEQC'],
                        DATABASE_FILE = DATABASE_FILE,
                        LOG_FILE = LOG_FILE,
                        START_FILE = startfilename)
            else:
                WQXtoPandas(
                    settingsDict['Input file'],
                    charDict,
                    outputPath = settingsDict['Path to output directory'],
                    outputDirName = settingsDict['Name of output directory'],
                    fromFile = True,
                    RUN_PHREEQC = RUN_PHREEQC,
                    bracket_charge_balance=bracket_charge_balance,
                    PHREEQC_PATH = settingsDict['Path to PHREEQC'],
                    DATABASE_FILE = DATABASE_FILE,
                    LOG_FILE = LOG_FILE,
                    START_FILE = startfilename)
        elif (settingsDict['Input method'] == '2'):
            #   We will use a list of sites from a NWIS XML file and query these
            #   sites for water quality data
            #First extract site list from XML file
            try:
                siteList = extractSitesFromXML(settingsDict['Input file'])
            except IOError:
                print("Problem extracting sites from XML file "+settingsDict['Input file']+" check to see if file name is correct and file is in right location.")
                return -1
            charList = []
            #collect list of characteristics to query
            for key in charDict.keys():
                charList.append(str(key))
            if len(siteList)>autosplitnum:
                #We have too long of a list and should split into multiple queries
                n_groups = int(ceil(len(siteList)/float(autosplitnum)))
                for i in range(n_groups):#this doesn't work for even division cases
                    shortList = siteList[i*autosplitnum:i*autosplitnum+autosplitnum]
                    queryText = querySiteList(shortList, charList)
                    if (queryText != None):
                        WQXtoPandas(
                            queryText,
                            charDict,
                            outputPath = settingsDict['Path to output directory'],
                            outputDirName = settingsDict['Name of output directory'],
                            fromFile = False,
                            RUN_PHREEQC = RUN_PHREEQC,
                            bracket_charge_balance=bracket_charge_balance,
                            PHREEQC_PATH = settingsDict['Path to PHREEQC'],
                            DATABASE_FILE = DATABASE_FILE,
                            splittag = '.'+str(i),
                            LOG_FILE = LOG_FILE,
                            START_FILE = startfilename)
            else:
                #get html for query
                queryText = querySiteList(siteList, charList)
                if (queryText != None):
                    WQXtoPandas(
                        queryText,
                        charDict,
                        outputPath = settingsDict['Path to output directory'],
                        outputDirName = settingsDict['Name of output directory'],
                        fromFile = False,
                        RUN_PHREEQC = RUN_PHREEQC,
                        bracket_charge_balance=bracket_charge_balance,
                        PHREEQC_PATH = settingsDict['Path to PHREEQC'],
                        DATABASE_FILE = DATABASE_FILE,
                        LOG_FILE = LOG_FILE,
                        START_FILE = startfilename)
        elif (settingsDict['Input method'] == '3'):
            #   We will use a list of sites from a text file and query these
            #   sites for water quality data
            #First extract site list from text file
            try:
                siteList = extractSitesFromText(settingsDict['Input file'])
            except IOError:
                print("Problem extracting sites from text file "+settingsDict['Input file']+" check to see if file name is correct and file is in right location.")
                return -1
            if (siteList != -1):
                charList = []
                #collect list of characteristics to query
                for key in charDict.keys():
                    charList.append(str(key))
                if len(siteList)>autosplitnum:
                    #We have too long of a list and should split into multiple queries
                    n_groups = int(ceil(len(siteList)/float(autosplitnum)))
                    for i in range(n_groups):
                        shortList = siteList[i*autosplitnum:i*autosplitnum+autosplitnum]
                        queryText = querySiteList(shortList, charList)
                        if (queryText != None):
                                     WQXtoPandas(
                                         queryText,
                                         charDict,
                                         outputPath = settingsDict['Path to output directory'],
                                         outputDirName = settingsDict['Name of output directory'],
                                         fromFile = False,
                                         RUN_PHREEQC = RUN_PHREEQC,
                                         bracket_charge_balance=bracket_charge_balance,
                                         PHREEQC_PATH = settingsDict['Path to PHREEQC'],
                                         splittag = '.'+str(i),
                                         DATABASE_FILE = DATABASE_FILE,
                                         LOG_FILE = LOG_FILE,
                                         START_FILE = startfilename)
                else:
                    #get html for query
                    queryText = querySiteList(siteList, charList)
                    if (queryText != None):
                        WQXtoPandas(
                            queryText,
                            charDict,
                            outputPath = settingsDict['Path to output directory'],
                            outputDirName = settingsDict['Name of output directory'],
                            fromFile = False,
                            RUN_PHREEQC = RUN_PHREEQC,
                            bracket_charge_balance=bracket_charge_balance,
                            PHREEQC_PATH = settingsDict['Path to PHREEQC'],
                            DATABASE_FILE = DATABASE_FILE,
                            LOG_FILE = LOG_FILE,
                            START_FILE = startfilename)
            else:
                print("Problem obtaining site list.")
        else:
            print(('Problem with "Input Method" of start file: ' + settingsDict['Input method']))

    except IOError:
        print("Problem reading start file.  Check file name.")


#Run as script
if __name__=="__main__":
    #pull in name of start file
    startfilename = sys.argv[1]
    if len(sys.argv)>2:
        autosplitnum = sys.argv[2]
        runWQXtoPandas(startfilename, autosplitnum=autosplitnum)
    else:
        runWQXtoPandas(startfilename)
