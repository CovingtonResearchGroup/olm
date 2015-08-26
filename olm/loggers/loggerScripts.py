#loggerScripts package
# main file contains functions that are useful in general for
# manipulation of data logger data

#from pandas import DateRange, depricated, use date_range
from pandas import datetools, DataFrame, Series, notnull
from matplotlib.dates import date2num

#accepts a list of logger DateFrame objects as first argument
def joinLoggers(loggerlist, how='inner', interpolate = False):
    #merge data from multiple loggers
    if type(loggerlist) == list:
        joined = loggerlist[0].join(loggerlist[1:], how=how)
        if interpolate:
            for col in joined.columns:
                filled_col = joined[col].interpolate()
                joined[col] = filled_col
        return joined
    else:
        print "Problem with input list: Need to input a list of DataFrame objects"
        return None

def joinAndResampleLoggers(loggerlist, interval, suffixes=[], how='inner', interpolate=False, limit=None):
    #If no suffixes were passed, create a list full of None values
    #  this keeps suffixes from being added in the code below
    if suffixes==[]:
        for i in arange(len(loggerlist)):
            suffixes.append(None)
    resampledList = []
    if type(loggerlist)==list:
        #print "Processing list type loggerlist..."
        for i,logger in enumerate(loggerlist):
            if suffixes[i]!=None:
                logger.columns+='_'+suffixes[i]
            resampledList.append(logger.resample(interval))
    elif type(loggerlist)==dict:
        #print "Processing dict type loggerlist..."
        for logger_key in loggerlist.keys():            
            logger = loggerlist[logger_key]
            if type(suffixes)==dict:
                if suffixes[logger_key]!=None:
                    logger.columns+='_'+suffixes[logger_key]
                resampledList.append(logger.resample(interval))
            else:
                print "Problem with suffixes. If loggerlist is a dict, suffixes also must be a dict."
                return None
    else:
        print "Problem with logger list: Need to input a list or dict of DataFrame or Series objects"
        return None
            
    for i, logger in enumerate(resampledList):
        if i==0:
            joined=logger
#            elif i==1:
#                joined=joined.join(logger, how=how, lsuffix='_'+suffixes[0], rsuffix='_'+suffixes[1])
#        elif i==3:
#            return joined
        else:
            joined=joined.join(logger, how=how)#, rsuffix='_'+suffixes[i])
    if interpolate:
        for col in joined.columns:
            print joined
            print col
            filled_col = joined[col].interpolate(limit=limit)
            joined[col] = filled_col
    return joined


def linear_correction(rawSeries, correctionSeries):
    #loop through correction series and calculate multiplying factors
    corrDict = {}
    for date, measurement in correctionSeries.iteritems():
        candidates = rawSeries.index[notnull(rawSeries)]
        index = candidates.searchsorted(date)
        if index > 0:
            asOfDate = candidates[index - 1]
            this_k = measurement/rawSeries[asOfDate]
            corrDict[asOfDate]= this_k
        else:
            asOfDate = rawSeries.index[0]
    if not rawSeries.index[0] in corrDict:
        corrDict[rawSeries.index[0]]=1
    if not rawSeries.index[-1] in corrDict:
        corrDict[rawSeries.index[-1]] = corrDict[asOfDate]
    k_series = Series(corrDict)
    joined_series = DataFrame({'raw':rawSeries,'k':k_series})
    joined_series.k = joined_series.k.interpolate()
    rawSeries = rawSeries*joined_series.k
    return rawSeries

#accepts a Series or DataFrame object as input
#def resampleHourly(ts):
#    start = ts.index[0]
#    start = start - datetime.timedelta(minutes=start.minute,
#                                       seconds=start.second)
#    end = ts.index[-1]
#    hourly = DateRange(start, end, offset=datetools.Hour())
#    grouped = ts.groupby(hourly.asof)
#    resampled = grouped.mean()
#    return resampled



#Function to extract logger data with same indicies as manual
#measurements both data sets are resampled on an hourly interval to
#assure alignment of indicies
# Arguments:
#    logger - a TimeSeries containing a single variable from the 
#             data logger along with time stamps
#    manual - a DataFrame or TimeSeries containing a single variable from the 
#             manual data set
#    value_name - a string to use as the base for the collumn labels.  
#                 Otherwise determined from input series.
#  Returns a DataFrame object containing both values and the aligned index
def manualCompare(logger, manual, value_name='', ltag='log', rtag='man'):
    if not(value_name==''):
        value_name += '_'
    logger = resampleHourly(logger)
    manual = resampleHourly(manual)
    wantidx = manual.index
    logger = logger[wantidx]
#    logger.name = value_name + '_logger'
#    manual.name = value_name + '_manual'
    joined = DataFrame({value_name+ltag:logger, value_name+rtag:manual})
#logger.align(manual, how='inner', lsuffix = '_log', rsuffix='_man')
 #   joined.rename(columns=[value_name+'_log', value_name+'_man' ])
    return joined
