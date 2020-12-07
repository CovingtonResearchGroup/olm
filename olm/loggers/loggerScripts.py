"""
Contains functions that are useful in general for manipulation of data logger data
"""

from pandas import DataFrame, Series, notnull, Timestamp
from matplotlib.dates import date2num
from numpy import arange

#accepts a list of logger DateFrame objects as first argument
def joinLoggers(loggerlist, how='inner', interpolate = False):
    """
    Joins the data from a list of logger DataFrame objects together.

    Parameters
    ----------
    loggerlist : list
        A list of logger DataFrame objects to be joined.
    how : string, optional
        How the two DataFrames are to be joined. Default is inner.
    interpolate : boolean, optional
        Determines whether empty rows are to be filled with data via interpolation. Uses Pandas Dataframe.interpolate(). Default = False
    
    Returns
    -------
    joined : pandas.DataFrame
      DataFrame of joined loggers.

    """
    #merge data from multiple loggers
    if type(loggerlist) == list:
        joined = loggerlist[0].join(loggerlist[1:], how=how)
        if interpolate:
            for col in joined.columns:
                filled_col = joined[col].interpolate()
                joined[col] = filled_col
        return joined
    else:
        print("Problem with input list: Need to input a list of DataFrame objects")
        return None

def joinAndResampleLoggers(loggerlist, interval, suffixes=[], how='inner', interpolate=False, limit=None):
    """
    Joins and resamples data from DataFrame objects provided in a list.

    Parameters
    ----------
    loggerlist : list
        List of logger pandas.core.dataframe.DataFrame objects to be joined.
    interval : string
        Pandas offset string (http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases) on which the DataFrames should be resampled (e.g. 'H'=hour, 'T'=minute, 'D'=day).
    suffixes : list
        A list of strings (same length as loggerlist) that contains suffixes to be applied to each logger. This is useful if multiple loggers have the same column names.
    how : string
        Method for joining loggers (default = 'inner').
    interpolate : boolean
        Whether data should be interpolated to fill gaps in rows (default=False).
    limit : int
        Maximum number of consecutive NaNs to fill if data are interpolated.

    Returns
    -------
    joined : pandas.core.dataframe.DataFrame
        DataFrame Object that contains joined DataFrames.

    """
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
            resampledList.append(logger.resample(interval).mean())
    elif type(loggerlist)==dict:
        #print "Processing dict type loggerlist..."
        for logger_key in list(loggerlist.keys()):            
            logger = loggerlist[logger_key]
            if type(suffixes)==dict:
                if suffixes[logger_key]!=None:
                    logger.columns+='_'+suffixes[logger_key]
                resampledList.append(logger.resample(interval).mean())
            else:
                print("Problem with suffixes. If loggerlist is a dict, suffixes also must be a dict.")
                return None
    else:
        print("Problem with logger list: Need to input a list or dict of DataFrame or Series objects")
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
#            print joined
#            print col
            filled_col = joined[col].interpolate(limit=limit)
            joined[col] = filled_col
    return joined


def linear_correction(rawSeries, correctionSeries):
    """
    Uses a Pandas Series of spot measured values to linearly correct time series data from a logger.

    Parameters
    ----------
    rawSeries : pandas.core.series.Series
        A Pandas Series that contains the time series data to be corrected.
    correctionSeries : pandas.core.series.Series
        A Pandas Series that contains spot measurement data that are to be used to correct rawSeries.

    Returns
    -------
    rawSeries : pandas.core.series.Series
        A corrected version of the rawSeries time series data.
    """
    #loop through correction series and calculate multiplying factors
    corrDict = {}
    for date, measurement in correctionSeries.items():
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

def manualCompare(logger, manual, value_name='', ltag='_log', mtag='_man'):
    """
    Function to extract logger data with same timestamps as manual measurements for comparison. Both data sets are resampled on an hourly interval to assure alignment of indicies.

    Parameters
    ----------
    logger : pandas.core.series.Series
        A Pandas TimeSeries containing a single column and time stamps as indices.
    manual : pandas.core.series.Series 
        A Pandas TimeSeries containing a single variable from the manual data set
    value_name : string 
        A string to use as the base for the collumn labels in the output DataFrame. 
    ltag : string
        A suffix to be added to the logger column name, or used as the logger column name if value_name is not set.
    mtag : string
        A suffix to be added to the manual measurement column name, or used as the manual measurement column name if value_name is not set.

    Returns
    -------
    joined : pandas.core.dataframe.DataFrame
        A DataFrame object containing values of manual measurements and corresponding values from the logger time series using the aligned index that is resampled to the hour.
    """
    if not(value_name==''):
        value_name += '_'
    logger = resampleHourly(logger)
    manual = resampleHourly(manual)
    wantidx = manual.index
    logger = logger[wantidx]
    joined = DataFrame({value_name+ltag:logger, value_name+rtag:manual})
    return joined


def shiftLogger(logger, shift_to, align_at_start = True):
    """
    Function to extract logger data with same timestamps as manual measurements for comparison. Both data sets are resampled on an hourly interval to assure alignment of indicies.

    Parameters
    ----------
    logger : pandas.core.series.Series or pandas.core.dataframe.Dataframe
        A Pandas TimeSeries or DataFrame containing time stamps as indices.
    shift_to : string 
        A string that contains the date and time that the logger series should be shifted to. By default this is the correct starting time (first time stamp) of the series.
    align_at_start : boolean
        If True, shift_to is assumed to represent the correct starting date for the series. If False, shift_to is assumed to represent the correct final date of the series. (default=True)

    Returns
    -------
    logger : pandas.core.series.Series or pandas.core.dataframe.DataFrame
        A Series or DataFrame object that contains the correct shifted time stamps.
    """
    bad_times = logger.index
    #align at starting time stamp
    if align_at_start:
        start_time = Timestamp(shift_to)
        dt = start_time - bad_times[0]
    #align at ending time stamp
    else:
        end_time = Timestamp(shift_to)
        dt = end_time - bad_times[-1]
    #shift index of original logger time series
    logger.index = logger.index + dt
    return logger
