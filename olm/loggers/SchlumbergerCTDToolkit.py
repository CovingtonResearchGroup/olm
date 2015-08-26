#Tools for reading and analysis of data from Schlumberger CTD Divers

from pandas import read_csv
from pandas import concat
from pandas import DataFrame

#read in the CSV file from a CTD diver and return a pandas DataFrame
def readCTD(csvfile):
    if csvfile.endswith('MON'):
        sep = '\s\s\s\s*'
    else:
        sep = ','
#    header = 0
    skiprows = 66 #this is somewhat weak, number of lines could change over time??
    index_col = 0
    names = ['Pressure', 'Temperature', 'Conductivity']
    parse_dates = True
    skip_footer = 1
    df = read_csv(csvfile, sep=sep, names=names, skiprows=skiprows, index_col=index_col, parse_dates=parse_dates, skip_footer=skip_footer)
    return df

#read in the CSV file from a CTD diver and return a pandas DataFrame
def readBaro(csvfile):
    if csvfile.endswith('MON'):
        sep = '\s\s\s\s*'
    else:
        sep = ','
#    header = 0
    skiprows = 54 #this is somewhat weak, number of lines could change over time??
    index_col = 0
    names = ['Pressure', 'Temperature']
    parse_dates = True
    skip_footer = 1
    df = read_csv(csvfile, sep=sep, names=names, skiprows=skiprows, index_col=index_col, parse_dates=parse_dates, skip_footer=skip_footer)
    return df


# concatCTD accepts a list of CTD DataFrames and concatenate them
#
# Arguments:
# dflist = a list of dataframes
#
# zero_shift = if set to True, the pressure values will be adjusted at the time
#              of each join, assuming that flow depth before and after the join
#              was equal.  If set to False, no adjustment will be made in      
#              pressure values
#
# n_to_average = number of data points to average before and after join in order
#                to determine data offset value for pressure
#
# offset_list = list of offsets to be applied manually to pressure data

def concatCTD(dflist, zero_shift = True, n_to_average = 5, offset_list=[], offset_dates = []):
    concatenated = None
    if zero_shift == False:
        #concatenate with no shifting
        #note: might want to add some capability to handle overlapping data
        concatenated = concat(dflist)
    else:
        if len(offset_list) > 0:
            #offset each data file by the value in offset list
            if len(offset_list) != len(dflist) - 1:
                print "Number of elements in offset_list must be one less than number of data files to concatenate"
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
