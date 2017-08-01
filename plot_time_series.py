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
from scipy.stats import t 
from atsr import *
from geolocation import getdata

#Paths;
figpath = '/group_workspaces/cems2/nceo_generic/CCI_LAND/figs/'
filepath = "/group_workspaces/cems/aerosol_cci/public/cci_products/AATSR_ORAC_v04-01/L3_MONTHLY"
outpath = "/home/users/fespir/Backup/CCI_LAND/pathinfo/"
#locate .nc files
list_of_paths = [] # empty list to write files ending in .nc to 


for root, dirs, files in os.walk(filepath): # locating files in the L3_Monthly directory ending in .nc 
    for file in files:
        if file.endswith("_seg.nc"):
            continue
        else:
            if file.endswith(".nc"):
                file_path = (os.path.join(root, file)) #gives the full root to the file path. e.g. /group_workspaces/cems/aerosol_cci/public/cci_products/AATSR_ORAC_v04-01/L3_MONTHLY/2008/200805-ESACCI-L3C_AEROSOL-AER_PRODUCTS-AATSR-ENVISAT-ORAC-MONTHLY-fv04.01.nc
                list_of_paths.append(file_path) # writes each full file path to the empty list defined above

#extracting the year and month from each file path name. E.g. 200805 2008 = year, 05 = month.          
list_of_times = []
list_of_months = [] 
list_of_years = []
list_of_files = []

#print atsr.aerosol_file_info(list_of_paths)


for item in list_of_paths:
    output = aerosol_file_info(item)
    print output
    list_of_months.append(output['MONTH'])
    list_of_years.append(output['YEAR'])
    list_of_files.append(output['filename'])
    ts = output['YEAR'] + (output['MONTH']-0.5)/12. #centre of the month
    list_of_times.append(ts)# writing each time to a new list


sortedtimes = np.argsort(list_of_times) #sorting list_of_times numerically
sortedpath = [list_of_paths[i] for i in sortedtimes] #ordering list_of_paths using corresponding sortedtimes array. 



AODs = [] # empty list where AOD's at -63.5 lon, -11.5 lat for each .nc file will be written to.
count1 = 0

# Extract info from .nc file
nan = [] 
counting = 0
for item in sortedpath:
    #print item
    values = getdata(item)
    print values 
    '''
    readfile = Dataset(item,mode='r') #read data 
    lons = readfile.variables['longitude'][:] 
    lats = readfile.variables['latitude'][:]
    mean_AOD_values = readfile.variables['AOD550_mean'][:,:]
    coordlon = np.where(lons == -63.5) # locate where in list of lons where value == -63.5 
    coordlat = np.where(lats == -11.5)
    AOD_at_point = mean_AOD_values[coordlat,coordlon]
    for i in AOD_at_point: #finding AOD values == nan
        if i < 0.01:
            nan.append(counting)#writing the indexes to a list 
    AOD_at_point = float(AOD_at_point) 
    AODs.append(AOD_at_point)
    counting = counting + 1
'''
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
#sys.exit("End")
#plotting monthly AOD mean
plt.xticks(month, xticks)
plt.plot(month, Averages, 'o-')
plt.title('Mean Monthly AOD 2002-2012')
plt.ylabel('AOD')
plt.xlabel('Month')
plt.show()
plt.savefig( figpath + 'monthly_AOD_mean_plot.pdf')

#plotting AOD anomalies


#plt.title('AOD anomalies 2002-2012')
#plt.xlabel('Year')
#plt.ylabel('AOD anomaly')
#A = np.vstack([X, np.ones(len(X))]).T
#m, c = np.linalg.lstsq(A, Y)[0]
#print m, c 
#bestfittxt =  '{:.5f}x + {:.10f}'.format(m, c)
#plt.plot(X, Y, 'o', label='Original data', markersize=5) 
#plt.plot(X, m*X + c, 'r', label='Fitted line ' + bestfittxt)

xs = np.array([1,5,6,7])
ys = np.array([1,4,6,6])
print xs,ys
As = np.vstack([xs, np.ones(len(xs))]).T
grad, inter = np.linalg.lstsq(As, ys)[0]
plt.plot(xs, ys, 'o', label='Original data', markersize=5) 
plt.plot(xs, grad*xs + inter, 'r', label='Fitted line ' + bestfittxt)
print grad, inter

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