from pandas import read_csv

def read_cat(cat_file):
    """
    Reads in data from a concatenated Cyclops file.

    Parameters
    ----------
    dat_file : string
        The name of the Cat.txt file to read.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing data from cat file.
    """
    skiprows = [0,1,2,3,4,6]
    index_col= 2
    parse_dates = True
    skipinitialspace =True
    df = read_csv(cat_file, skiprows=skiprows, index_col=index_col, parse_dates=parse_dates, skipinitialspace=skipinitialspace)
    return df

def read_single(txt_file):
    """
    Reads in data from a single non-concatenated Cyclops file.

    Parameters
    ----------
    txt_file : string
        The name of the Cat.txt file to read.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing data from cat file.
    """
    skiprows = 2
    index_col= 0
    parse_dates = True
    df = read_csv(txt_file, skiprows=skiprows, index_col=index_col, parse_dates=parse_dates)
    return df
