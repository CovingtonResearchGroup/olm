#Tools for reading and analysis of data from Schlumberger CTD Divers

from pandas import read_csv
from pandas import concat
from pandas import DataFrame

"""
Functions to read Schlumberger diver logger files.
"""

#read in the CSV file from a CTD diver and return a pandas DataFrame
def readCTD(csvfile):
    """
    Reads data from a CSV or MON file exported from a Schlumberger CTD Diver.

    Parameters
    ----------
    csv_file : string
        A string containing the file name of the CSV or MON file to be read.
       
    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing data from HOBO csv file.
    """
    if csvfile.endswith('MON'):
        sep = '\s\s\s\s*'
    else:
        sep = ','
#    header = 0
    skiprows = 66 #this is somewhat weak, number of lines could change over time??
    index_col = 0
    names = ['Pressure', 'Temperature', 'Conductivity']
    parse_dates = True
    skipfooter = 1
    df = read_csv(csvfile, sep=sep, names=names, skiprows=skiprows, index_col=index_col, parse_dates=parse_dates, skipfooter=skipfooter)
    return df

#read in the CSV file from a CTD diver and return a pandas DataFrame
def readBaro(csvfile):
    """
    Reads data from a CSV or MON file from a Schlumberger Baro Diver.

    Parameters
    ----------
    csv_file : string
        A string containing the file name of the CSV or MON file to be read.
       
    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing data from HOBO csv file.
    """
    if csvfile.endswith('MON'):
        sep = '\s\s\s\s*'
    else:
        sep = ','
#    header = 0
    skiprows = 54 #this is somewhat weak, number of lines could change over time??
    index_col = 0
    names = ['Pressure', 'Temperature']
    parse_dates = True
    skipfooter = 1
    df = read_csv(csvfile, sep=sep, names=names, skiprows=skiprows, index_col=index_col, parse_dates=parse_dates, skipfooter=skipfooter)
    return df


# 
#
# Arguments:
# dflist = a list of dataframes
#
# zero_shift = 
#
# #
# 

def concatCTD(dflist, zero_shift = True, n_to_average = 5, offset_list=[], offset_dates = []):
    """
    Accepts a list of CTD DataFrames and concatenates them.

    Parameters
    ----------
    dflist : list
        List of pandas.DataFrames to concatenate.
    zero_shift : boolean
        If set to True, the pressure values will be adjusted at the time of each join, assuming that flow depth before and after the join was equal.  If set to False, no adjustment will be made in pressure values. This is useful when downloading the logger may have resulted in a slightly different position in the water column. (Default = True)
    n_to_average : int
        Number of data points to average before and after join in order to determine data offset value for pressure
    offset_list : list
        List of offsets to be applied manually to pressure data.
    offset_dates : list
        List of datetime strings corresponding to manual offsets.

    Returns
    -------
    (concatenated : pandas.DataFrame, offset_list : pandas.DataFrame)
        A tuple is returned with the first item being a DataFrame object containing the concatenated data and the second item in the tuple being a DataFrame object containing offsets with datetimes of the offsets as an index.
        
    """
    concatenated = None
    if zero_shift == False:
        #concatenate with no shifting
        #note: might want to add some capability to handle overlapping data
        concatenated = concat(dflist)
    else:
        if len(offset_list) > 0:
            #offset each data file by the value in offset list
            if len(offset_list) != len(dflist) - 1:
                print("Number of elements in offset_list must be one less than number of data files to concatenate")
                return None
            else:
                for i, df in enumerate(dflist):
                    if i != 0: #skip first data frame
                        df['Pressure'] = df['Pressure'] + offset_list[i-1]
        else:
            for i, df in enumerate(dflist):
                if i != 0: #skip first data frame
                    #in tail/head we throw out last/first data point
                    #get average value from tail of previous data
                    tail_values = dflist[i-1]['Pressure'][-n_to_average-1:-1]
                    tail_average = tail_values.mean()
                    #get average value from head of following data
                    head_values = df['Pressure'][1:n_to_average+1]
                    head_average = head_values.mean()
                    delta = tail_average - head_average
                    offset_dates.append(df.index[0])
                    offset_list.append(delta)
                    df['Pressure'] = df['Pressure'] + delta
        concatenated = concat(dflist)
    offsets = DataFrame(offset_list, index=offset_dates)
    return (concatenated, offsets)
