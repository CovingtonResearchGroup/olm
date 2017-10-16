#Functions to handle reading and manipulating data from the WTW-pH data loggers

from pandas import read_csv
from pandas import concat
from pandas import DataFrame
from dateutil.parser import parse

def readpH(csvfile):
    """
    Reads data from a CSV file produced by WTW-pH data loggers.

    Parameters
    ----------
    csvfile : string
        Name of the CSV file to be read.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame object containing pH data.
    """
    sep = ';'
    names = ['SensorName', 'SerNo', 'DateTime', 'pH', 'Status-ph', 'Temperature', 'Status-temp', 'MeasID', 'CalStatus', 'CalProtID', 'blank']
    index_col = 3
    skiprows = 5
    parse_dates = True
    #function to convert commas to decimal points and strings to floats in data
    decimal_conversion = lambda x: float(x.replace(",","."))
    converters = {'pH':decimal_conversion, 'Temperature':decimal_conversion}
    temp_df = read_csv(csvfile, sep=sep, skiprows=skiprows, names=names,index_col=index_col, date_parser=euroParser, converters=converters)
    df = temp_df[['Temperature', 'pH']]
    return df

def euroParser(datestring):
    """
    Date parser for European style dates separated by '.' characters.

    Parameters
    ----------
    datestring : string
        String containing Euro style date

    Returns
    -------
    dateTime : datetime object
    """
    splitdate = datestring.split('.')
    if len(splitdate) == 3:
        newdatestring = splitdate[1] +'.' + splitdate[0] +'.' + splitdate[2]
    else:
        newdatestring = datestring
    dateTime = parse(newdatestring)
    return dateTime

