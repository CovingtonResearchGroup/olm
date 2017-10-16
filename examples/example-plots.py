#Once you have run the example input file
#'runWQXtoPandas-Sample-Input-File.xls' then you should be able to
#make some initial plots of the data using the code below.  This
#provides only a short example of what can be done.

import pickle as pickle
from pylab import *

#Read in the Pandas Panel that contains the data for one site from the pickled output file.
site_pnl = pickle.load(open('alaska/USGS-15024800/USGS-15024800-Panel.pkl', 'rb'))
#The panel contains multiple layers (think of them as spreadsheet sheets). One layer contains the data, the others contain meta data.  You could possibly use the meta data to filter what data you actually want to analyze.
#For now, we will just grab the data (raw #'s)
data = site_pnl['data']
#Read in PHREEQC output
data_phrq = pickle.load(open('alaska/USGS-15024800/USGS-15024800-PHREEQC.pkl', 'rb'))
#Join two dataframes and add a suffix to any duplicate column names from PHREEQC
data = data.join(data_phrq, rsuffix='_phrq')

#These are the column labels we have available
print("Available data:")
print((data.columns.values))
print(" ")
#Grab temperature data
temp = data['Temperature, water']
#Plot it out as a time series
figure()
temp.plot()
ylabel("Water temperature (C)")
show()
#Get alkalinity
alk = data['Alkalinity, total']
#Get mean discharges
discharge = data['Stream flow, mean. daily']
#Plot alkalinity vs. discharge on a loglog plot with squares as symbols
figure()
loglog(discharge, alk, 's')
xlabel('Discharge (cfs)')
ylabel('Alkalinity (as mg/L CaCO3)')
show()
#We can also do some simple stats
print("Alkalinity stats:")
print((alk.describe()))
print(" ")

#Plot SI Calcite vs. Discharge
figure()
SI = data['SI_Calcite']
semilogx(discharge, SI, 'o')
xlabel('Discharge (cfs)')
ylabel('SI Calcite')
show()

#What about site data
site_data = pickle.load(open('alaska/USGS-15024800/USGS-15024800-Site-Description.pkl', 'rb'))
area = site_data['drain_area_va']
timezone = site_data['tz_cd']
name = site_data['station_nm']
lat = site_data['dec_lat_va']
longitude = site_data['dec_long_va']

print('Station name = ', name)
print('Latitude = ', lat)
print('Longitude = ', longitude)
print('Time Zone = ', timezone)
print('Drainage basin area = ', area)
