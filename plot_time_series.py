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

#Paths;
figpath = '/group_workspaces/cems2/nceo_generic/CCI_LAND/figs/'
path = "/group_workspaces/cems/aerosol_cci/public/cci_products/AATSR_ORAC_v04-01/L3_MONTHLY/"
outpath = "/group_workspaces/cems2/nceo_generic/CCI_LAND/output_data"

#Check if save file containing these files already exists
aerfile = outpath+'aodfiles.json'


if not os.path.isfile(aerfile):
    #All files are located in extract_files that are inputs to fetch_aerosol_file_attributes
    #that outputs a structure containing all file times.
    aerosol_files = fetch_aerosol_file_attributes(extract_files(path,"04.01.nc")  )
    # Writing JSON data
    print aerosol_files
    #sys.exit()
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

file_AOD_lon_lat = outpath + 'aodlonlat.txt'
print file_AOD_lon_lat

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
    nan = [] 
    counting = 0
    for item in sortedpath:
        readfile = Dataset(item,mode='r')
        lons = readfile.variables['longitude'][:] 
        lats = readfile.variables['latitude'][:]
        mean_AOD_values = readfile.variables['AOD550_mean'][:,:]
        coordlon = np.where(lons == -63.5) # locate where in list of lons where value == -63.5 
        coordlat = np.where(lats == -11.5)
        AOD_at_point = mean_AOD_values[coordlat,coordlon]
        print AOD_at_point   
        #np.where(data_reading < 0.01):  
        for i in AOD_at_point: #finding AOD values == nan
            if i < 0.01:
                nan.append(counting)#writing the indexes to a list 
        AOD_at_point = float(AOD_at_point) 
        AODs.append(AOD_at_point)
        counting = counting + 1
    with open(file_AOD_lon_lat, 'w') as writing:
        for value in AODs:       
            writing.write("%s\n" % value )
      
with open(file_AOD_lon_lat) as reading:
    content = reading.readlines()
    print content
    
#reading = open(file_AOD_lon_lat.txt, 'r')
#print reading
sys.exit()

#print AODs


ordered_times = [list_of_times[i] for i in sortedtimes]#odering list_of_times numericaly  
#removing the nan values from the list of times and list of AODs
x = [i for j, i in enumerate(ordered_times) if j not in nan]
y = [i for j, i in enumerate(AODs) if j not in nan]

#sys.exit("stopping")

#plot the AODs for each month from 2002-2012
plt.plot(x, y, 'o-')
plt.title('Monthly AOD 2002-2012')
plt.ylabel('AOD')
plt.xlabel('Month')
plt.show()
plt.savefig( figpath + 'monthly_AOD_2002_2012.pdf')

AOD_values = np.array(y) #converting AODs to an array for the def function below
X = np.array(x)

print X
#sys.exit("exiting")
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
#removing dates when AOD == nan: 
months = [i for j, i in enumerate(months_ordered) if j not in nan]
for item in months:
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


monthly_retrievals = []
monthly_retrievals.extend((len(Jan), len(Feb),len(Mar),len(Apr),len(May),len(Jun),len(Jul),len(Aug),len(Sep),len(Oct),len(Nov),len(Dec)))
monthly_retrievlas = np.array(monthly_retrievals) 

#compute anomalie for each AOD value by substracting the monthly mean
anomfinal = [] 
for i in range(len(AOD_values)):
    anom2 = AOD_values[i] - Averages[int(months[i])-1]
    anomfinal.append(anom2)
Y = np.array(anomfinal)
print 'anomalies', Y
#plotting AOD monthly average
xticks = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'] # my x ticks 
month = linspace(1,12,12) # list from 1 to 12 

plt.xticks(month,xticks)
plt.plot(month,monthly_retrievals, 'o')
plt.title('Monthly Retrievals')
plt.show()


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
plt.xlabel('Year')
plt.ylabel('AOD anomaly')
A = np.vstack([X, np.ones(len(X))]).T
m, c = np.linalg.lstsq(A, Y)[0]
y_x0= m*X[0] + c
bestfittxt =  '{:.5f}x + {:.10f}'.format(m, y_x0)
plt.plot(X, Y, 'o', label='Original data', markersize=5) 
plt.plot(X, m*X + c, 'r', label='Fitted line ' + bestfittxt)
#intercept at X[0]:

sys.exit("end")

mu = sum(Y)/len(Y)
ttest = stats.ttest_1samp(Y, mu)
print 't-statistic = %6.3f pvalue = %6.4f' % ttest #equal_var = False)
print stats.t.ppf(1-0.05, len(Y))

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
