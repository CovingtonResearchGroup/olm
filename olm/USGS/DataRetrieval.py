#Function to download USGS Daily Average Discharge values given date and the location text
"""
Contains functions used to download data from USGS databases.
"""

from lxml import etree
import lxml.html, requests
try:
    from urllib.parse import quote #could eventually rework to use only requests
except ImportError:
    from urllib import quote
from io import StringIO
from pandas import read_csv, DataFrame, to_datetime
#import requests
#import os


def GetSiteData(location):
    """
    Retrieves meta data about a Water Quality portal site using the site identifier.


    Parameters
    ----------
    location : string
        Full site number.

    Returns
    -------
    siteDF : pandas.DataFrame
        Returns a pandas DataFrame object that contains all of the site meta data.

    Notes
    -----

    """
    BASEURL = 'https://www.waterqualitydata.us/data/Station/search?siteid='
    queryURL = BASEURL + location + '&mimeType=csv&Zip=no'
    #Need to skip header, which is hopefully uniform across USGS queries
    siteDF = read_csv(queryURL, sep=',')
    siteDF = siteDF.iloc[0]#change axis so that we only have key-data pairs
    return siteDF


def GetNWISSiteData(location):
    """
    Retrieves meta data about a USGS site using the full site identifier.


    Parameters
    ----------
    location : string
        Full USGS site number starting with 'USGS-' or the bare integer number of a USGS site.

    Returns
    -------
    siteDF : pandas.DataFrame
        Returns a pandas DataFrame object that contains all of the site meta data from an expanded USGS site data query.  Data is indexed using field labels given in USGS file (see Notes).

    Notes
    -----
    Pair of keys and descriptions from USGS site meta data.

    agency_cd       -- Agency

    site_no         -- Site identification number

    station_nm      -- Site name

    site_tp_cd      -- Site type

    lat_va          -- DMS latitude

    long_va         -- DMS longitude

    dec_lat_va      -- Decimal latitude

    dec_long_va     -- Decimal longitude

    coord_meth_cd   -- Latitude-longitude method

    coord_acy_cd    -- Latitude-longitude accuracy

    coord_datum_cd  -- Latitude-longitude datum

    dec_coord_datum_cd -- Decimal Latitude-longitude datum

    district_cd     -- District code

    state_cd        -- State code

    county_cd       -- County code

    country_cd      -- Country code

    land_net_ds     -- Land net location description

    map_nm          -- Name of location map

    map_scale_fc    -- Scale of location map

    alt_va          -- Altitude of Gage/land surface

    alt_meth_cd     -- Method altitude determined

    alt_acy_va      -- Altitude accuracy

    alt_datum_cd    -- Altitude datum

    huc_cd          -- Hydrologic unit code

    basin_cd        -- Drainage basin code

    topo_cd         -- Topographic setting code

    instruments_cd  -- Flags for instruments at site

    construction_dt -- Date of first construction

    inventory_dt    -- Date site established or inventoried

    drain_area_va   -- Drainage area

    contrib_drain_area_va -- Contributing drainage area

    tz_cd           -- Time Zone abbreviation

    local_time_fg   -- Site honors Daylight Savings Time

    reliability_cd  -- Data reliability code

    gw_file_cd      -- Data-other GW files

    nat_aqfr_cd     -- National aquifer code

    aqfr_cd         -- Local aquifer code

    aqfr_type_cd    -- Local aquifer type code

    well_depth_va   -- Well depth

    hole_depth_va   -- Hole depth

    depth_src_cd    -- Source of depth data

    project_no      -- Project number

    """
    if (location[:5] == 'USGS-'):
        sitenum = location[5:]
    else:
        sitenum = location
    BASEURL = 'https://waterservices.usgs.gov/nwis/site/?site='
    queryURL = BASEURL + sitenum + '&siteOutput=expanded'
    #Need to skip header, which is hopefully uniform across USGS queries
    skiprows = list(range(0,59))
    skiprows.append(60)
    siteDF = read_csv(queryURL, sep='\t', skiprows=skiprows)
    siteDF = siteDF.iloc[0]#change axis so that we only have key-data pairs
    return siteDF

def querySiteList(siteList, charList):
    BASE_URL = 'https://www.waterqualitydata.us/Result/search?'
    queryText = BASE_URL + 'siteid='
    #add sites to query
    for site in siteList:
        #check for USGS prefixes (are there others?  EPA?)
        ##if not(site.startswith('USGS-')):
        ##    site = 'USGS-' + site
        #add this site to list with trailing semi-colon
        queryText += site + ';'
    #remove final semi-colon
    queryText = queryText[:-1]
    #add characteristics to query
    queryText += '&characteristicName='
    for characteristic in charList:
        queryText += characteristic + ';'
    #remove trailing semi-colon
    queryText = queryText[:-1]
    #add mime type
    queryText += '&mimeType=xml'
    #convert query string to url special characters
    queryText = quote(queryText, safe="/&=:?")
    return queryText


def GetDailyDischarge(location, date):
    """
    Retrieve daily average discharge value from USGS database for given date and USGS site.

    Parameters
    ----------
    location : string
        Full USGS site number starting with 'USGS-' or a string that just contains the bare integer number of a USGS site.
    date : string
        String containing the date for which discharge will be retrieved.  Should be given as YYYY-MM-DD.

    Returns
    -------
    data : dict {'discharge':float, 'quality':string, 'name':string}
        Returns a dicionary that contains three items, the average discharge value for the site and date given, the quality code assigned to that discharge value, and the name of the site.

    Notes
    -----
    Currently hard-wired to retrieve USGS pcode 00060, daily discharge in cfs.

    """
    #construct url for discharge query
    BASE_URL = 'https://waterservices.usgs.gov/nwis/dv?format=waterml,1.1'
    #query discharge and read into xml parser
    #pull site number out of location text
    #Check to see if location contains 'USGS-' or is just the bare number
    if (location[:5] == 'USGS-'):
        site_number = location[5:]
    else:
        site_number = location
    #construct html address for query
    query_html = BASE_URL + '&sites=' + site_number + '&startDT='+date+'&endDT='+date
    #read in xml file through html query
    try:
        print("Discharge query html: ",query_html)
        r = requests.get(query_html)
        #qtree = etree.parse(r.raw)
        root = etree.fromstring(r.content)
    except IOError:
        print("Problem retrieving discharge value (IOError).")
        return -1
    #parse xml file to pull out discharge and quality code
    #root = qtree.getroot()
    #get namespace map
    NSMAP = root.nsmap
    NS1 = "{%s}" % NSMAP['ns1']
    tsString = "timeSeries[@name='USGS:"+site_number+":00060:00003']"
    ts = root.find(NS1+tsString)
    if (ts == None):
        #there is no time series data for this site and date
        return None
    sourceInfo = ts.find(NS1+"sourceInfo")
    name = sourceInfo.findtext(NS1+"siteName")
    values = ts.find(NS1+"values")
    value = values.find(NS1+"value")
    if (value == None):
        #there is no discharge data for this site and date
        return None
    q = value.text
    quality_code = value.get("qualifiers")
    #return discharge and quality code
    data = {'discharge':q, 'quality':quality_code, 'name':name}
    return data


def GetDailyDischargeRecord(location, start_date, end_date=None):
    """
    Retrieve daily average discharge values from USGS database for given date range and USGS site.

    Parameters
    ----------
    location : str
        Full USGS site number starting with 'USGS-' or a string that just contains the bare integer number of a USGS site.
    start_date : str
        String containing the beginning date in the range for which discharge will be retrieved.  Should be given as YYYY-MM-DD.
    end_date : str (optional)
        String containing the ending date in the range for which discharge will be retrieved.  Should be given as YYYY-MM-DD.  If not provided then data will be retrieved up to the current date.

    Returns
    -------
    data : pandas dataframe
        Returns a Pandas dataframe with an index of the date, a column 'discharge' of discharge values, and a column 'quality' of the USGS quality rating.

    Notes
    -----
    Currently hard-wired to retrieve USGS pcode 00060, daily discharge in cfs.

   """
    #construct url for discharge query
    BASE_URL = 'https://waterservices.usgs.gov/nwis/dv?format=waterml,1.1'
    #query discharge and read into xml parser
    #pull site number out of location text
    #Check to see if location contains 'USGS-' or is just the bare number
    if (location[:5] == 'USGS-'):
        site_number = location[5:]
    else:
        site_number = location
    #construct html address for query
    if end_date==None:
        query_html = BASE_URL + '&sites=' + site_number + '&startDT='+start_date
    else:
        query_html = BASE_URL + '&sites=' + site_number + '&startDT='+start_date+'&endDT='+end_date
    #read in xml file through html query
    try:
        r = requests.get(query_html)
        #qtree = etree.parse(r.raw)
        root = etree.fromstring(r.content)
    except IOError:
        print("Problem retrieving discharge value (IOError).")
        return -1
    #parse xml file to pull out discharge and quality code
    #root = qtree.getroot()
    #get namespace map
    NSMAP = root.nsmap
    NS1 = "{%s}" % NSMAP['ns1']
    tsString = "timeSeries[@name='USGS:"+site_number+":00060:00003']"
    ts = root.find(NS1+tsString)
    if (ts == None):
        #there is no time series data for this site and date
        return None
    sourceInfo = ts.find(NS1+"sourceInfo")
    name = sourceInfo.findtext(NS1+"siteName")
    values = ts.find(NS1+"values")
    value_list = values.findall(NS1+"value")
    if (values == None):
        #there is no discharge data for this site and date
        return None
    q=[]
    quality_code=[]
    date_list = []
    for value in value_list:
        q.append(float(value.text))
        quality_code.append(value.get("qualifiers"))
        date_list.append(value.get("dateTime")[:10])
    #4/24/14 ended coding here.  need to write into dataframe
    data = DataFrame({'discharge':q, 'quality':quality_code}, index=to_datetime(date_list))
    #return discharge and quality code data frame
    return data
