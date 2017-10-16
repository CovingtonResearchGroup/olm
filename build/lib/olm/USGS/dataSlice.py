##################################################################
## functions to slice desired data from data dicts produced by  ##
## loadWaterQualityData.py functions                            ##
##################################################################

# dataDict = a general data dict from loadWaterQualityData
# iterKeyList = the keys needed to get to the level of the dict where we 
#               want to iterate and find data (e.g. the site level).  If
#               set to none, it is assumed that the dict is passed at 
#               iteration level
# dataKeyList = an ordered list of the keys (one or more) needed 
#               to get from the iteration level to the desired data

def extractValues(dataDict, dataKeyList=None, iterKeyList=None):
    if (iterKeyList == None):
        iterDict = dataDict
    else:
        tempIterDict = dataDict
        for key in iterKeyList:
            tempIterDict = tempIterDict[key]
        iterDict = tempIterDict
    values = []
    gaps = []
    for subDictKey in iterDict.keys():
        subDict = iterDict[subDictKey]
        temp = subDict
        isGap = False
        for key in dataKeyList:
            if (key in temp):
                temp = temp[key]
            else:
                #this is a data gap
                gaps.append(subDictKey)
                isGap = True
                break
        if not(isGap):
            values.append(temp)
    return {'values':values, 'gaps':gaps}

def extractValueSet(dataDict, dataKey2DList, iterKeyList = None, getKeys = False, getGaps = False): 
    if (iterKeyList == None):
        iterDict = dataDict
    else:
        tempIterDict = dataDict
        for key in iterKeyList:
            tempIterDict = tempIterDict[key]
        iterDict = tempIterDict
    valuesList = []
    gapsList = []
    #create blank 2D lists with the proper number of top level elements
    usedIterKeys = []
    for list in enumerate(dataKey2DList):
        valuesList.append([])
        gapsList.append([])
    for subDictKey in iterDict.keys():
        subDict = iterDict[subDictKey]
        isGap = False
        variableList = []
        for variableIdx, dataKeyList in enumerate(dataKey2DList):
            temp = subDict
            for key in dataKeyList:
                if (key in temp):
                    temp = temp[key]
                else:
                    #this is a data gap
                    gapsList[variableIdx].append(subDictKey)
                    isGap = True
                    break
            if not(isGap):
                variableList.append(temp)     
        if not(isGap):            
            usedIterKeys.append(subDictKey)
            for i, variable in enumerate(variableList):
                valuesList[i].append(variable)
    setDict = {}
    if (getKeys):
        # we want the keys
        setDict['keys'] = usedIterKeys
    if (getGaps):
        # we want the gaps
        setDict['gaps'] = gapsList
    if ( len(setDict) > 0 ):
        #either keys or gaps were added so return as dict with desired components
        setDict['values'] = valuesList
        return setDict
    else:
        # we only want the 2D list of values
        return valuesList
