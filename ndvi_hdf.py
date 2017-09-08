#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 13:35:31 2017

@author: fespir
"""


import numpy as np
from file_search import *
from pyhdf.SD import SD, SDC
import pprint
import sys
from atsr import *
from geo_ndvi import *
from doy_to_month import *
import json
import csv 
from geolocation import *
import matplotlib.pyplot as plt 
from math import *

figpath = "/group_workspaces/cems2/nceo_generic/CCI_LAND/figs/"
path = "/group_workspaces/cems2/nceo_generic/satellite_data/modis_c6/ndvi/"
outpath = "/group_workspaces/cems2/nceo_generic/CCI_LAND/output_data/"

NDVI_out = outpath + "NDVI_files.json"

F = (extract_files(path,'.hdf'))

try:    
    file = open(NDVI_out,'r')
except IOError: # writing to file if it doesn't exist
    fileinfo = fetch_NDVI_file_attributes(F)
    file = open(NDVI_out, 'wb')
    with file as f:
        json.dump(fileinfo, f)  

    
with open(NDVI_out,'r') as reading:
    filedata = json.load(reading)
    Year = filedata['Year']
    Month = filedata['Month']
    Day = filedata['Day']
    Time = filedata['Times']

M = []
T = []
Y = []
D = []
for i in range(len(Month)):
    M.append(int(Month[i]))
    Y.append(int(Year[i]))
    T.append(float(Time[i]))
    D.append(int(Day[i]))
 
SORT = np.argsort(T)
YEAR = np.asarray(Y)[SORT]
DAY = np.asarray(D)[SORT]
TIME = np.asarray(T)[SORT]
MONTH = np.asarray(M)[SORT]
FILE = np.asarray(F)[SORT]



lat_bnds, lon_bnds = [-11, -9], [-65.5, -63.5]
#lat_bnds, lon_bnds = [-12, -11], [-64, -63]
mean_ndvi_out = outpath + "mean_ndvi_out.csv"
'''
try:
    file = open(mean_ndvi_out, "r")
except IOError:     
    ndvi_means = []
    for item in F:
        readdata = get_NDVI_readings(item,lat_bnds,lon_bnds,'CMG 0.05 Deg Monthly NDVI',10000., 0.)
        mean_ndvi = grid_average(readdata)
        ndvi_means.append(mean_ndvi)
    with open(mean_ndvi_out, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in ndvi_means:
            writer.writerow([val]) 
'''            
ndvi_means = []
for item in F:
        readdata = get_NDVI_readings(item,lat_bnds,lon_bnds,'CMG 0.05 Deg Monthly NDVI',10000., 0.)
        mean_ndvi = grid_average(readdata)
        ndvi_means.append(mean_ndvi)
with open(mean_ndvi_out, "w") as output:
    writer = csv.writer(output, lineterminator='\n')
    for val in ndvi_means:
        writer.writerow([val])
ndvi_vals = []
with open(mean_ndvi_out,"r") as r:
        reader = csv.reader(r,delimiter=' ')
        for row in reader:   
            values = float((",".join(row)))
            ndvi_vals.append(values)
            
NDVIs = np.asarray(ndvi_vals)
print NDVIs

#calculating monthly average and writing to list
data = monthly_average(MONTH,NDVIs)
monthly_retrievals = data[1]
NDVI_average = data[0]

anomalies = [] 
#calculating anomalies for each of the precipitation values
for i in range(len(MONTH)):
    anom = NDVIs[i] - NDVI_average[int(MONTH[i])-1]
    anomalies.append(anom)

print anomalies
#converting plotting list to array
anomalies = np.asarray(anomalies)
monthly_retrievals = np.asarray(monthly_retrievals)
Y1 = np.array(anomalies)

xtick = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
months = linspace(1,12,12)
year_list_NDVI = linspace(2000,2017,18)


def scatter_LOBF(Time, Data, xs, ylab, xlab, scatterlabel, figtitle,figname):
#plotting AOD anomalies
    fig = plt.figure(figsize=(16,4))
    ax = fig.add_subplot(111)
    plt.title(figtitle, fontsize = 14)
    plt.xlabel(xlab,fontsize = 14 )
    plt.ylabel(ylab, fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xticks(xs,fontsize = 14)
    plot = plt.plot(Time, Data, 'o', label=scatterlabel, markersize=5) 
    matrix = np.vstack([Time, np.ones(len(Time))]).T
    m, c = np.linalg.lstsq(matrix, Data)[0]

    #finding intercept at the time of the first AOD reading
    y_x0= m*Time[0] + c

    #.5f =  5 significant figures
    bestfittxt =  '{:.5f}x + {:.10f}'.format(m, y_x0)
    
    #plotting line of best fit 
    bestfit = plt.plot(Time, m*Time + c, 'r', label='Fitted line ' + bestfittxt)

    #calculating the mean, variance, T statistic and p value
    m, v, s, k = stats.t.stats(5, moments='mvsk')
    n, (smin, smax), sm, sv, ss, sk = stats.describe(Data)
    sstr = 'mean = %6.4f, variance = %6.4f'
    tt = (sm-m)/np.sqrt(sv/float(n))  # t-statistic for mean
    pval = stats.t.sf(np.abs(tt), n-1)*2  # two-sided pvalue = Prob(abs(t)>tt)
    print 't-statistic = %6.3f pvalue = %6.4f' % (tt, pval)
    print 'distribution: ', sstr %(m, v)
    print 'sample: ', sstr %(sm, sv)
    plt.legend()
    plt.show()
    plt.savefig( figpath + figname)
    return plot, bestfit 


def scatter(tick_spacing,xticks, retrievals, ylab, xlab, figtitle, figname):
    #number of retrievlas for each month from 2002-2012
    plt.xticks(tick_spacing,xticks, fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.plot(tick_spacing,retrievals, 'o-')
    plt.title(figtitle, fontsize = 14)
    plt.ylabel(ylab, fontsize = 14)
    plt.xlabel(xlab, fontsize = 14)
    plot = plt.show()
    plt.savefig(figpath+figname)
    return plot

monthly_ave_plot = scatter(months, xtick, NDVI_average, 'monthly NDVI average', 'Year', 'Monthly NDVI average Rondonia','ndvi_monthly_average_largeA.png')
Anomaly_plot = scatter_LOBF(TIME, Y1, year_list_NDVI,'NDVI anomalies','Time', 'NDVI','NDVI Anomalies 2000-2017 - Rondonia','anomalies_NDVI_kar.png')
#Procedure to obtain latitude and longitude values for
#MODIS NDVI data


'''
        

fileinfo = fetch_NDVI_file_attributes(extract_files(path,'.hdf'))
with open(NDVI_out, "w") as output:
        writer = csv.writer(output)
        for key, value in fileinfo.items():
            writer.writerow([key, value])

with open(NDVI_out,"r") as r:
        reader = csv.reader(r, quoting=csv.QUOTE_NONNUMERIC)
        reading = dict(reader)
        Day = reading["Day"]
        Month = reading["Month"]
        Year = reading["Year"]
        Times = reading["Times"]
        File = reading["File"] 
print type(Times)
'''
