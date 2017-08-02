"""
Created on Wed Jul 26 12:08:40 2017

@author: fespir
"""

"""
NAME;
Monthly time series of AOD over the 2002-2012 period for the Rondonia region 

PURPOSE;
This program processes monthly AOD data at 550nm for the Rondonia region (11S,63W) from 2002 - 2012.
Monthly time series plots are produced (AOD vs TIME) where TIME = YEAR + MONTH/12. 

DESCRIPTION;
The processor uses information from 114 monthly AATSR datasets from 2002 - 2012.

The programme saves each .nc file to a lsit, extracts the date from the file path name and orders the file paths
by the date. Next, the .nc files are read and the AOD, lat and lon variables are extracted. The mean AOD for the 
Rondonia region is then plotted 114 times for each date on a AOD vs TIME plot. 

A second plot of AOD vs MONTH is produced by calculating the mean AOD for each month over the 10 year period.

Finally, the mean AOD for each month is taken from the 114 AOD values to determine the anomaly for each point and
plotted as AOD vs TIME. 

INPUT;
NetCDF files containing the AOD550 variables for each pixel location

OUTPUT;
Time series plot 

AUTHOR;
Freya Espir
"""
import json
import pandas as pd
from numpy import * 
from math import *
import os,sys
from netCDF4 import Dataset
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
import fnmatch
import matplotlib.pyplot as plt 
import itertools
from operator import itemgetter 
from scipy import stats
from itertools import repeat
from scipy.stats import t 
from atsr import *
from file_search import *
from geolocation import *
import csv

#Paths;
figpath = '/group_workspaces/cems2/nceo_generic/CCI_LAND/figs/'
path = "/group_workspaces/cems/aerosol_cci/public/cci_products/AATSR_ORAC_v04-01/L3_MONTHLY/"
outpath = "/group_workspaces/cems2/nceo_generic/CCI_LAND/output_data"

#Check if save file containing these files already exists
aerfile = outpath+'aodfiles.json'


if os.path.isfile(aerfile):
    #All files are located in extract_files that are inputs to fetch_aerosol_file_attributes
    #that outputs a structure containing all file times.
    aerosol_files = fetch_aerosol_file_attributes(extract_files(path,"04.01.nc")  )
    # Writing JSON data
    #print aerosol_files
    with open(aerfile, 'w') as f:
        json.dump(aerosol_files, f)
#sys.exit()
# Reading data back
with open(aerfile, 'r') as f:
    aerosol_files = json.load(f)
    list_of_times  = aerosol_files['time']
    list_of_months = aerosol_files['month']
    list_of_years  = aerosol_files['year']
    list_of_files  = aerosol_files['file']




#sys.exit('stop')
sortedtimes = np.argsort(list_of_times) #sorting list_of_times numerically
sortedpath = [list_of_files[i] for i in sortedtimes] #ordering list_of_paths using corresponding sortedtimes array. 

#print sortedpath

file_AOD_lon_lat = outpath + 'aodlonlat.json'

#if os.path.isfile(file_AOD_lon_lat):
    #gettingdata = getdata(item,'longitude', 'latitude', 'AOD550_mean',-63.5,-11.5)
    #print gettingdata['AATSR_READING']
    #sys.exit()
    #listofaods.append(gettingdata)#['AATSR_READING'])
    #json2pd = pd.Series(gettingdata).to_json(orient='values')
    #with open(file_AOD_lon_lat, 'w') as f:
    # json.dump(json2pd,f) 
    #listofaods = []
    #aerosol_readings = json.load(opening)
    #data_reading = aerosol_readings['AATSR_READING']

if not os.path.isfile(file_AOD_lon_lat):
  
    AODs = [] 
    for item in sortedpath:
        readfile = Dataset(item,mode='r')
        lons = readfile.variables['longitude'][:] 
        lats = readfile.variables['latitude'][:]
        mean_AOD_values = readfile.variables['AOD550_mean'][:,:]
        coordlon = np.where(lons == -63.5) # locate where in list of lons where value == -63.5 
        coordlat = np.where(lats == -11.5)
        AOD_at_point = mean_AOD_values[coordlat,coordlon]
        AOD_at_point = float(AOD_at_point) 
        AODs.append(AOD_at_point)  
    with open (file_AOD_lon_lat, 'w') as outfile:
        json.dump(AODs,outfile)

    #nans = np.argwhere(np.isnan(AODs)) 
    #removednan = [i for j, i in enumerate(AODs) if j not in nans]
    #np.asarray(removednan)
    #np.savetxt(file_AOD_lon_lat, removednan, delimiter="'")   
with open(file_AOD_lon_lat) as reading:
    content = json.load(reading)
    #print content
    
#converting lists to arrays
times = np.asarray(list_of_times)
months = np.asarray(list_of_months)
years = np.asarray(list_of_years)
files = np.asarray(list_of_files)
aod = np.asarray(AODs)


#index of sorted times 
SORT = np.argsort(times)

#sorting lists usings SORT index
TIMES = times[SORT]
AOD = aod[SORT]
MONTHS = months[SORT]
YEARS = years[SORT]
FILES = files[SORT]



#removing nan values
AOD_no_nan = outpath + 'AOD_no_nan.json' 

#
#print filtered 


nonan = ~np.isnan(AOD)
y = AOD[nonan]
x = TIMES[nonan]
MONTHS = MONTHS[nonan]
YEARS = YEARS[nonan]
FILES = FILES[nonan]


y = array2json(y)
x = array2json(x)
MONTHS = array2json(MONTHS)
YEARS = array2json(YEARS)
FILES = array2json(FILES)

#filtered = remove_nan(AOD,TIMES,MONTHS,YEARS,FILES)
#y = pd.Series(y).to_json(orient='values')

with open(AOD_no_nan, 'w') as towrite:
    json.dump(y, towrite)
    
with open(AOD_no_nan, 'r') as toread:
    data = json.load(toread)
    print data
sys.exit()

#Monthly AOD average
count = 1
AOD_average = []
monthly_retrievals = []
month_list = linspace(1,12,12)

for i in range(len(month_list)):
    location = np.where(MONTHS == count)
    location = np.asarray(location)
    AODMONTH = y[location]
    total = sum(AODMONTH)
    length = AODMONTH.size
    average = total/length
    AOD_average.append(average)
    monthly_retrievals.append(length)
    count = count+1 



       
    

#compute anomalie for each AOD value by substracting the monthly average
anomalies = [] 
for i in range(len(MONTHS)):
    anom = y[i] - AOD_average[int(MONTHS[i])-1]
    anomalies.append(anom)

#converting plotting list to array
anomalies = np.asarray(anomalies)
monthly_retrievals = np.asarray(monthly_retrievals)


#plotting AOD monthly average
xticks = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'] # my x ticks 


#plot the AODs for each month from 2002-2012
plt.plot(x, y, 'o-')
plt.title('Monthly AOD 2002-2012')
plt.ylabel('AOD')
plt.xlabel('Month')
plt.show()
plt.savefig( figpath + 'monthly_AOD_2002_2012.pdf')


plt.xticks(month_list,xticks)
plt.plot(month_list,monthly_retrievals, 'o')
plt.title('Monthly Retrievals')
plt.show()


#plotting monthly AOD mean
plt.xticks(month_list, xticks)
plt.plot(month_list, AOD_average, 'o-')
plt.title('Mean Monthly AOD 2002-2012')
plt.ylabel('AOD')
plt.xlabel('Month')
plt.show()
plt.savefig( figpath + 'monthly_AOD_mean_plot.pdf')

#plotting AOD anomalies


plt.title('AOD anomalies 2002-2012')
plt.xlabel('Year')
plt.ylabel('AOD anomaly')
A = np.vstack([x, np.ones(len(x))]).T
m, c = np.linalg.lstsq(A, anomalies)[0]
y_x0= m*x[0] + c
print y_x0
bestfittxt =  '{:.5f}x + {:.10f}'.format(m, y_x0)
plt.plot(x, anomalies, 'o', label='Original data', markersize=5) 
plt.plot(x, m*x + c, 'r', label='Fitted line ' + bestfittxt)
#intercept at X[0]:


mu = sum(anomalies)/len(anomalies)
ttest = stats.ttest_1samp(anomalies, mu)
print 't-statistic = %6.3f pvalue = %6.4f' % ttest #equal_var = False)
print stats.t.ppf(1-0.05, len(anomalies))

m, v, s, k = stats.t.stats(5, moments='mvsk')
n, (smin, smax), sm, sv, ss, sk = stats.describe(Y)
sstr = 'mean = %6.4f, variance = %6.4f'
print 'distribution: ', sstr %(m, v)
print 'sample: ', sstr %(sm, sv)
tt = (sm-m)/np.sqrt(sv/float(n))  # t-statistic for mean
pval = stats.t.sf(np.abs(tt), n-1)*2  # two-sided pvalue = Prob(abs(t)>tt)
print 't-statistic = %6.3f pvalue = %6.4f' % (tt, pval)
plt.legend()
plt.show()
plt.savefig( figpath + 'AOD_anomalies.pdf')
