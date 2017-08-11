#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 14:10:39 2017
NAME;
Monthly time series of precipitation from May 2000 - March 2017 for the Rondonia region 

PURPOSE;
This program processes monthly precipitation data (mm/hr) for the Rondonia region (11S,63W).
Plots:
1. Precipitation rate vs month
2. Number of precipitation retrievals over the 17 year period vs month
3. Mean precipitation per month vs month
4. Monthly precipitation anomalies vs year 


DESCRIPTION;
The processor uses information from 203 monthly TRMM datasets from 2000 - 2017.

1. The programme saves each .nc file to a list, extracts the date from the file path name and orders the file paths
by the date. 

2. .nc files are read and the precipitation rate, lat and lon variables for eaach .nc file are extracted. 

3. Calculates mean precipitation per month over the time period

4. calculates anomaly for each precipitation measurement using the monthly mean


INPUT;
NetCDF files containing the precipitation variables for each pixel location

OUTPUT;
Four time series plots 

AUTHOR;
Freya Espir

"""

import json
import pandas as pd
from numpy import * 
from math import *
import sys
import netCDF4 
from netCDF4 import Dataset
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt 
from operator import itemgetter 
from scipy import stats
from itertools import repeat
from scipy.stats import t 
from atsr import *
from file_search import *
from geolocation import *
import csv
import StringIO

#Paths;
figpath = "/group_workspaces/cems2/nceo_generic/CCI_LAND/figs/"
path = "/group_workspaces/cems2/nceo_generic/satellite_data/precipitation/TRMM/monthly"
outpath = "/group_workspaces/cems2/nceo_generic/CCI_LAND/output_data/"

#Check if save file containing these files already exists
precip_outpath = outpath+'precipitation.txt'

try:    
    file = open(precip_outpath,'r')
except IOError: # writing to file if it doesn't exist
    precipitation = fetch_aerosol_file_attributes(extract_files(path,"_TRMM_3B43.7.nc"))
    file = open(precip_outpath, 'w')
    with file as f:
        json.dump(precipitation, f) 

# Reading data back
with open(precip_outpath, 'r') as f:
    attributes = json.load(f)
    dates  = attributes['time']
    months = attributes['month']
    years  = attributes['year']
    files  = attributes['file']


#list to array
date = np.asarray(dates)
m = np.asarray(months)
yr = np.asarray(years)
f = np.asarray(files)


#sorting files numerically by the date
SORT = np.argsort(date)
TIMES = date[SORT]
FILES = f[SORT]
YEARS = yr[SORT]
MONTHS = m[SORT] 
          
mean_precip_out = outpath = 'precip_means.csv'

#latitude and longitude cooridnates for Rondonia region 
lat_bnds, lon_bnds = [-12, -11], [-64, -63]


try: 
    file = open(mean_precip_out, 'r')
except IOError:
    all_precip_means = []
    for item in FILES:
        precip = getdata_coordinategrid(item, 'nlon','nlat',lat_bnds,lon_bnds,'precipitation')                       
        mean_precip = grid_average(precip)
        all_precip_means.append(mean_precip)
    with open(mean_precip_out, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in all_precip_means:
            writer.writerow([val]) 

#ndvi = outpath + 'ndvi_mean.csv'

precip_vals = []
with open(mean_precip_out,"r") as r:
        reader = csv.reader(r,delimiter=' ')
        for row in reader:   
            values = float((",".join(row)))
            precip_vals.append(values)
precipitations = np.asarray(precip_vals)


#calculating monthly average and writing to list
data = monthly_average(MONTHS,precipitations)
monthly_retrievals = data[1]
precip_average = data[0]

anomalies = [] 

#calculating anomalies for each of the precipitation values
for i in range(len(MONTHS)):
    anom = precipitations[i] - precip_average[int(MONTHS[i])-1]
    anomalies.append(anom)

#print anomalies

#converting plotting list to array
anomalies = np.asarray(anomalies)
monthly_retrievals = np.asarray(monthly_retrievals)

Y1 = np.array(anomalies)

#stadard deviation for anomalies 
print 'anomalies standard dev', np.std(Y1, axis = 0)

# my x ticks 
xticks = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'] 

year_list = linspace(2000,2018,19)
month_list = linspace(1,12,12) 


#plot the precipitation for each month from 2002-2012
fig = plt.figure(figsize=(20,4))
ax = fig.add_subplot(111)
plt.yticks(fontsize = 14)
plt.xticks(year_list, fontsize = 12)
ax.plot(TIMES,precipitations, 'o-')
plt.title('Monthly Precipitation May 2000 - March 2017 \n', fontsize = 16)
plt.ylabel('Mean Precipitation(mm/hr)', fontsize = 14)
plt.xlabel('Month', fontsize = 14)
#plt.xlim([200000,200250])
plt.show()
#plt.savefig( figpath + 'monthly_AOD_2002_2012.pdf')

#number of retrievlas for each month from 2002-2012
plt.xticks(month_list,xticks, fontsize = 14)
plt.yticks(fontsize = 14)
plt.plot(month_list,monthly_retrievals, 'o')
plt.title('Monthly Retrievals \n', fontsize = 16)
plt.ylabel('Total Number of Retrievals', fontsize = 14)
plt.xlabel('Month', fontsize = 14)

plt.show()


#plotting monthly pmean
plt.xticks(month_list, xticks, fontsize = 14)
plt.yticks(fontsize = 14)
plt.plot(month_list, precip_average, 'o-')
plt.title('Average Precipitation per Month, May 2000- March 2017 \n', fontsize = 16)
plt.ylabel('Average Precipitation (mm/hr)', fontsize = 14)
plt.xlabel('Month', fontsize = 14)
plt.show()
#plt.savefig( figpath + 'monthly_AOD_mean_plot.pdf')


#plotting precipitation anomalies
plt.title('Precipitation Anomalies May 2000 - March 2017 \n', fontsize = 16)
plt.xlabel('Year',fontsize = 14 )
plt.ylabel('Precipitation Anomaly (mm/hr)', fontsize = 14)
plt.yticks(fontsize = 14)
plt.xticks(fontsize = 10)
matrix = np.vstack([TIMES, np.ones(len(TIMES))]).T
m, c = np.linalg.lstsq(matrix, Y1)[0]

#finding intercept at the time of the first AOD reading
y_x0= m*TIMES[0] + c

#.5f =  5 significant figures
bestfittxt =  '{:.5f}x + {:.10f}'.format(m, y_x0)

#plotting the points
plt.plot(TIMES, Y1, 'o', label='Original data', markersize=5) 

#plotting line of best fit 
plt.plot(TIMES, m*TIMES + c, 'r', label='Fitted line ' + bestfittxt)

m, v, s, k = stats.t.stats(5, moments='mvsk')
n, (smin, smax), sm, sv, ss, sk = stats.describe(precipitations)
sstr = 'mean = %6.4f, variance = %6.4f'
print 'distribution: ', sstr %(m, v)
print 'sample: ', sstr %(sm, sv)
tt = (sm-m)/np.sqrt(sv/float(n))  # t-statistic for mean
pval = stats.t.sf(np.abs(tt), n-1)*2  # two-sided pvalue = Prob(abs(t)>tt)
print 't-statistic = %6.3f pvalue = %6.4f' % (tt, pval)
       