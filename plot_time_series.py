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
outpath = "/group_workspaces/cems2/nceo_generic/CCI_LAND/output_data/"

#Check if save file containing these files already exists
aerfile = outpath+'aodfiles.json'

try: 
    file = open(aerfile, 'r')
except IOError:
    file = open(aerfile, 'w')
    aerosol_files = fetch_aerosol_file_attributes(extract_files(path,"04.01.nc"))
    with file as f:
        json.dump(aerosol_files, f)


# Reading data back
with open(aerfile, 'r') as f:
    aerosol_files = json.load(f)
    list_of_times  = aerosol_files['time']
    list_of_months = aerosol_files['month']
    list_of_years  = aerosol_files['year']
    list_of_files  = aerosol_files['file']



#converting lists to arrays
times = np.asarray(list_of_times)
months = np.asarray(list_of_months)
years = np.asarray(list_of_years)
files = np.asarray(list_of_files)


#sorting lists usings SORT index
SORT = np.argsort(list_of_times)
TIMES = times[SORT]
MONTHS = months[SORT]
YEARS = years[SORT]
FILES = files[SORT]

    
file_AOD_lon_lat = outpath + 'aodlonlat.csv'


#extracting AOD measurement, latitude and longitude information for a given coordinate 
#writing to a file if it does not already exist 

try: 
    file = open(file_AOD_lon_lat, 'r')
except IOError:
    aodmeans = []
    for item in FILES:
        gettingdata = getdata(item,'longitude', 'latitude', 'AOD550_mean',-63.5,-11.5)
        aodmeans.append(gettingdata)# getdata(filepath, longitude name in .nc file, latitude name in .nc file, variable name in .nc file, lon coordinate, lat coordinate)
    with open(file_AOD_lon_lat, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in aodmeans:
            writer.writerow([val]) 
             

       
AODs = []
with open(file_AOD_lon_lat,"r") as r:
        reader = csv.reader(r,delimiter=' ')
        for row in reader:   
            AODvals = float((",".join(row)))
            AODs.append(AODvals)
#reading file and extracting data 
  

#converting AODs to array
AOD_VALUE = np.asarray(AODs)

#removing nan values
#~np.isnan finds index of all non nan values in AODs list
nonan = np.where(~np.isnan(AODs))

#removing item from list if the corresponding AOD value is = nan
A = AOD_VALUE[nonan]
S = SORT[nonan]
T = TIMES[nonan]
M = MONTHS[nonan]
Y = YEARS[nonan]
F = FILES[nonan]

#converting A array to json format
Y = array2json(A)

AOD_no_nan = outpath + 'AOD_no_nan.json' 

try: 
    file = open(AOD_no_nan, 'r')
except IOError:
    file = open(AOD_no_nan, 'w')
    with file as towrite:
        json.dump(Y, towrite)
    
with open(AOD_no_nan, 'r') as toread:
    data = json.load(toread)

    
month_list = linspace(1,12,12)

#Calculating monthly AOD average
count = 1
AOD_average = []
monthly_retrievals = []

for i in range(len(month_list)):
    location = np.where(M == count)
    location = np.asarray(location)
    AODMONTH = A[location]
    total = sum(AODMONTH)
    length = AODMONTH.size
    average = total/length
    AOD_average.append(average)
    monthly_retrievals.append(length)
    count = count+1 


#computing anomalie for each AOD value by substracting the monthly average
anomalies = [] 

for i in range(len(M)):
    anom = A[i] - AOD_average[int(M[i])-1]
    anomalies.append(anom)


#converting plotting list to array
anomalies = np.asarray(anomalies)
monthly_retrievals = np.asarray(monthly_retrievals)

Y1 = np.array(anomalies)
print 'standard dev', np.std(Y1, axis = 0)
# my x ticks 
xticks = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'] 
year_list = linspace(2002,2012,11)

#plot the AODs for each month from 2002-2012
plt.plot(T, A, 'o-')
plt.title('Monthly AOD 2002-2012', fontsize = 14)
plt.yticks(fontsize = 14)
plt.xticks(fontsize = 14)
plt.ylabel('AOD', fontsize = 14)
plt.xlabel('Month', fontsize = 14)
plt.show()
plt.savefig( figpath + 'monthly_AOD_2002_2012.pdf')

#number of retrievlas for each month from 2002-2012
plt.xticks(month_list,xticks, fontsize = 14)
plt.yticks(fontsize = 14)
plt.plot(month_list,monthly_retrievals, 'o')
plt.title('Monthly Retrievals', fontsize = 14)
plt.ylabel('Total Number of Retrievals', fontsize = 14)
plt.xlabel('Month', fontsize = 14)

plt.show()


#plotting monthly AOD mean
plt.xticks(month_list, xticks, fontsize = 14)
plt.yticks(fontsize = 14)
plt.plot(month_list, AOD_average, 'o-')
plt.title('Mean Monthly AOD 2002-2012', fontsize = 14)
plt.ylabel('AOD', fontsize = 14)
plt.xlabel('Month', fontsize = 14)
plt.show()
plt.savefig( figpath + 'monthly_AOD_mean_plot.pdf')

#plotting AOD anomalies
fig = plt.figure(figsize=(16,4))
ax = fig.add_subplot(111)
plt.title('AOD anomalies 2002-2012', fontsize = 14)
plt.xlabel('Year',fontsize = 14 )
plt.ylabel('AOD Anomaly', fontsize = 14)
plt.yticks(fontsize = 14)
plt.xticks(year_list,fontsize = 14)
matrix = np.vstack([T, np.ones(len(T))]).T
m, c = np.linalg.lstsq(matrix, Y1)[0]

#finding intercept at the time of the first AOD reading
y_x0= m*T[0] + c

#.5f =  5 significant figures
bestfittxt =  '{:.5f}x + {:.10f}'.format(m, y_x0)

#plotting the points
plt.plot(T, Y1, 'o', label='Original data', markersize=5) 

#plotting line of best fit 
plt.plot(T, m*T + c, 'r', label='Fitted line ' + bestfittxt)

#calculating the mean, variance, T statistic and p value

m, v, s, k = stats.t.stats(5, moments='mvsk')
n, (smin, smax), sm, sv, ss, sk = stats.describe(A)
sstr = 'mean = %6.4f, variance = %6.4f'
print 'distribution: ', sstr %(m, v)
print 'sample: ', sstr %(sm, sv)
tt = (sm-m)/np.sqrt(sv/float(n))  # t-statistic for mean
pval = stats.t.sf(np.abs(tt), n-1)*2  # two-sided pvalue = Prob(abs(t)>tt)
print 't-statistic = %6.3f pvalue = %6.4f' % (tt, pval)
plt.legend()
plt.show()
plt.savefig( figpath + 'AOD_anomalies.pdf')
