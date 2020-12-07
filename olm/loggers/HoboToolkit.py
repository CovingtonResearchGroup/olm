from pandas import read_csv,DataFrame



def read_hobo_csv(csv_file, all_columns=False):
    """
    Reads data from a csv file exported from HOBOware.

    Parameters
    ----------
    csv_file : string
        A string containing the file name of the csv file to be read.
    all_columns : boolean (optional)
        Determines whether to read in all columns or just ones that we search for and relabel (RH, DewPt, Abs Pres, Temp, Attached, Stopped, Connected, EOF, Cond High Range, Cond Low Range, DO). Default = False
       
    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing data from HOBO csv file.
    """
    skiprows=1
    index_col=1
    parse_dates=True
    df = read_csv(csv_file, skiprows=skiprows, index_col=index_col, parse_dates=parse_dates)
    #Convert column names into something nicer
    columns = df.columns
    rename_dict = {}
    cond_count =0
    for label in columns:
        #By default, use old name
        new_name=label
        if all_columns==False:
            wantcol=False
        else:
            wantcol=True
        if label=='#':
            if all_columns:
                new_name = 'Num'
        if 'RH' in label:
            new_name = 'RH'
            wantcol=True
        if 'DewPt' in label:
            new_name = 'DewPt'
            wantcol = True
        if 'Abs Pres' in label:
            new_name = 'Pressure'
            wantcol=True
        if 'Temp' in label:
            new_name = 'Temp'
            wantcol=True
        if 'Detached' in label:
            if all_columns:
                new_name = 'Detached'
        if 'Attached' in label:
            if all_columns:
                new_name = 'Attached'
        if 'Stopped' in label:
            if all_columns:
                new_name = 'Stopped'
        if 'Connected' in label:
            if all_columns:
                new_name = 'Connected'
        if 'End Of File' in label:
            if all_columns:
                new_name = 'EOF'
        if 'Low Range' in label:
            new_name = 'CondLow'
            cond_count+=1
            wantcol=True
        if 'High Range' in label:
            new_name = 'CondHigh'
            cond_count+=1
            wantcol=True
        if 'Full Range' in label:
            new_name = 'CondFull'
            cond_count+=1
            wantcol=True
        if 'DO conc' in label:
            new_name = 'DO'
            wantcol=True
        if wantcol==True:
            rename_dict[label]=new_name
    #If there is only one conductivity column, we'll label it as 'Cond'
    if cond_count==1:
        old_names = list(rename_dict.keys())
        for old_name,new_name in rename_dict.items():
            if 'Cond' in new_name:
                cond_key = old_name
        rename_dict[cond_key] = 'Cond'
    df = df.rename(columns=rename_dict)
    if not(all_columns):
        #Trim out unwanted columns
        s_dict = {}
        for col in rename_dict.values():
            s = df[col]
            s_dict[col] = s
        df = DataFrame(s_dict)
    return df
