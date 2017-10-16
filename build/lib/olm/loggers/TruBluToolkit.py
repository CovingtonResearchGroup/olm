#Tools for reading and analysis of data from TruBlu data loggers 

from pandas import read_csv
from pandas import concat
from pandas import DataFrame
import os

"""
Functions to read TruBlu logger files.
"""

#read in the CSV file from a TruBlu logger and return a pandas DataFrame
def readTruBlu(csvfile):
    """
    Reads data from a CSV file exported from a TruBlu logger.

    Parameters
    ----------
    csv_file : string
        A string containing the file name of the CSV or MON file to be read.
       
    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing data from a TruBlu csv file.
    """
    sep = ','
    header = 0
    skiprows = 16 #this is somewhat weak, number of lines could change over time??
	# Definitely weak.  Probably an automated read to csv header would be better
    index_col = 3
    #names = ['ID','Name','Address','Time of Acquisition','Elapsed(Sec)','Level(PSI)','Temperature (\'C)','Battery Voltage(Volt)','Supply Voltage(Volt)','Scan No','blank']
    parse_dates = True
    #skip_footer = 1
    #print(csvfile)
    #df = read_csv(csvfile, sep=sep, names=names, skiprows=skiprows, index_col=index_col, parse_dates=parse_dates)
    
    try:
        if os.stat(csvfile).st_size > 0:
            df = read_csv(csvfile, sep=sep, skiprows=skiprows, header=header, index_col=index_col, parse_dates=parse_dates)
            return df
        else:
            print((csvfile + " is empty"))
    except OSError:
        print((csvfile + " does not exist"))


