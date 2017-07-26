
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
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


#Paths;
figpath = '/group_workspaces/cems2/nceo_generic/CCI_LAND/figs/'
filepath = "/group_workspaces/cems/aerosol_cci/public/cci_products/AATSR_ORAC_v04-01/L3_MONTHLY"

#locate .nc files
list_of_paths = [] # empty list to write files ending in .nc to 
for root, dirs, files in os.walk(filepath): # locating files in the L3_Monthly directory ending in .nc 
    for file in files:
        if file.endswith(".nc"):
            file_path = (os.path.join(root, file)) #gives the full root to the file path. e.g. /group_workspaces/cems/aerosol_cci/public/cci_products/AATSR_ORAC_v04-01/L3_MONTHLY/2008/200805-ESACCI-L3C_AEROSOL-AER_PRODUCTS-AATSR-ENVISAT-ORAC-MONTHLY-fv04.01.nc
            list_of_paths.append(file_path) # writes each full file path to the empty list defined above


#extracting the year and month from each file path name. E.g. 200805 2008 = year, 05 = month.          
list_of_times = []
list_of_months = [] 
list_of_years = []
list_of_files = []

for item in list_of_paths:
    getfilepath = item.split('/') # seperating the file path by '/' and writing each section as a new element in a list. E.g.   #print getfilepath
    # /AATSR_ORAC_v04-01/L3_MONTHLY/2008/200805-ESACCI-L3C_AEROSOL-AER_PRODUCTS-AATSR-ENVISAT-ORAC-MONTHLY-fv04.01.nc   
    # becomes a new list [AATSR_ORAC_v04-01,L3_MONTHLY,2008,200805-ESACCI-L3C_AEROSOL-AER_PRODUCTS-AATSR-ENVISAT-ORAC-MONTHLY-fv04.01.nc]
    # Date and month info is found in the last item in this list
    filename = getfilepath[(len(getfilepath)-1)] # last element in list
    #Extract the year
    YYYY = filename[:4]
    YEAR = float(YYYY)
    #Extract the month
    MM = filename[4:6]
    #print MM
    MONTH = float(MM)
    #Time series type plot
    list_of_months.append(MONTH)
    list_of_years.append(YEAR)
    list_of_files.append(filename)
    ts = YEAR + (MONTH-0.5)/12. #centre of the month
    list_of_times.append(ts)# writing each time to a new list

sortedtimes = np.argsort(list_of_times) #sorting list_of_times numerically
sortedpath = [list_of_paths[i] for i in sortedtimes] #ordering list_of_paths using corresponding sortedtimes array. 
x = [list_of_times[i] for i in sortedtimes]
X = np.array(x)
print X
print 'paths', len(list_of_paths)


AODs = [] # empty list where AOD's at -63.5, -11.5 for each .nc file will be written to.

# Extract info from .nc file
for item in sortedpath:
    readfile = Dataset(item,mode='r') #read data 
    lons = readfile.variables['longitude'][:] 
    lats = readfile.variables['latitude'][:]
    mean_AOD_values = readfile.variables['AOD550_mean'][:,:]
    coordlon = np.where(lons == -63.5) # locate where in list of lons where value == -63.5 
    coordlat = np.where(lats == -11.5)
    AOD_at_point = mean_AOD_values[coordlon,coordlat] # retrieve AOD value at this coordinate
    AOD_at_point = float(AOD_at_point) 
    AODs.append(AOD_at_point)


# plot the AODs for each month from 2002-2012
plt.plot(x, AODs, 'o-')
plt.title('Monthly AOD 2002-2012')
plt.ylabel('AOD')
plt.xlabel('Month')
plt.show()
plt.savefig( figpath + 'monthly_AOD_2002_2012.pdf')

AOD_values = np.array(AODs) #converting AODs to an array for the def function below

#Monthly AOD average
Jan = []
Feb = []
Mar = []
Apr = []
May = []
Jun = []
Jul = []
Aug = []
Sep = []
Oct = []
Nov = []
Dec = []
count = 0  
#List of months is in 1,2,3,4 etc format. Where 1 = Jan, 2 = Feb etc. This for loop identifies where 
#e.g. list_of_months ==1 and writes that index to a new list, in this case Jan. (could use np.where()?)

#ordering list_of_months using corresponding sortedtimes array. 
months_ordered = [list_of_months[i] for i in sortedtimes]

for item in months_ordered:
    if item == 1: # January 
        Jan.append(count) # e.g. 1st, 3rd and 5th element in list_of _months =1, so Jan = [0,2,4]
    if item == 2:
        Feb.append(count)
    if item == 3:
        Mar.append(count)
    if item == 4:
        Apr.append(count)
    if item == 5:
        May.append(count)
    if item == 6:
        Jun.append(count)
    if item == 7:
        Jul.append(count)
    if item == 8:
        Aug.append(count)
    if item == 9:
        Sep.append(count)
    if item == 10:
        Oct.append(count)
    if item == 11:
        Nov.append(count)
    if item == 12:
        Dec.append(count)
    count = count+1

Averages = []
#seperate AOD values into different months
def computeaverage(month,AOD):
    AOD_month = itemgetter(*month)(AOD) #Extracting AOD value for each month from list of all AOD values
    average = sum(AOD_month)/len(AOD_month) # avearge AOD for each month
    Averages.append(average) # writing monthly averages to a list 
    return AOD_month #return AOD values for each month 

#functions to be read into def comuteaverage(month,AOD) which computes monthly average and writes this to the Averages list.
m1 = computeaverage(Jan,AOD_values) # Jan = index
m2 = computeaverage(Feb, AOD_values)
m3 = computeaverage(Mar, AOD_values)
m4 = computeaverage(Apr, AOD_values)
m5 = computeaverage(May, AOD_values)
m6 = computeaverage(Jun, AOD_values)
m7 = computeaverage(Jul, AOD_values)
m8 = computeaverage(Aug, AOD_values)
m9 = computeaverage(Sep, AOD_values)
m10 = computeaverage(Oct, AOD_values)
m11 = computeaverage(Nov, AOD_values)
m12 = computeaverage(Dec, AOD_values)


#compute anomalie for each AOD value by substracting the monthly mean
anomfinal = [] 
for i in range(len(AOD_values)):
    anom2 = AOD_values[i] - Averages[int(months_ordered[i])-1]
    anomfinal.append(anom2)
Y = np.array(anomfinal)
#plotting AOD monthly average
xticks = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'] # my x ticks 
month = linspace(1,12,12) # list from 1 to 12 
 
#plotting monthly AOD mean
plt.xticks(month, xticks)
plt.plot(month, Averages, 'o-')
plt.title('Mean Monthly AOD 2002-2012')
plt.ylabel('AOD')
plt.xlabel('Month')
plt.show()
plt.savefig( figpath + 'monthly_AOD_mean_plot.pdf')

#plotting AOD anomalies
plt.title('AOD anomalies 2002-2012')
plt.xlabel('Month')
plt.ylabel('AOD anomaly')
plt.plot(x, anomfinal, marker = 'o', color = 'black',
        markerfacecolor='blue', markersize=6)    

#calculate linear least squares fit 
A = np.vstack([X, np.ones(len(X))]).T
m, c = np.linalg.lstsq(A, Y)[0]
print (m,c)
bestfittxt =  '{:.5f}x {:.5f}'.format(m, c)
#plot line and label inc. equation
plt.plot(X, m*X + c, 'r', label= 'Fitted line '+ bestfittxt)
plt.legend()
plt.show()
plt.savefig( figpath + 'AOD_anomalies.pdf')
