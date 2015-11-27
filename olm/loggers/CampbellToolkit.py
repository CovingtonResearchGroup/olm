from pandas import read_csv

def read_dat(dat_file):
    """
    Reads in data from a .dat file produced by Campbell's data logging software.

    Parameters
    ----------
    dat_file : string
        The name of the .dat file to read.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing data from .dat file.
    """
    skiprows = [0,2,3]
    index_col= 0
    parse_dates = True
    df = read_csv(dat_file, skiprows=skiprows, index_col=index_col, parse_dates=parse_dates, na_values='NAN')
    return df
