#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
NAME;
Script to process NDVI data for a coodinate grid and calculate anomalies, number of ndvi readings for each month and the average ndvi reading for each month across the data time period.

PURPOSE;
This program processes monthly NDVI data for the Rondonia region  (9-11S,63.5-65.5W).
Outputs lists for the follwing plots:
1. NDVI vs month
2. Number of NDVI retrievals over the 17 year period vs month
3. Average monthly NDVI vs month
4. Monthly precipitation anomalies vs year 


DESCRIPTION;
The processor uses information from 209 monthly MODIS datasets from 2000 - 2017.

1. The programme saves each .hdf file to a list, extracts the year and day of year from the file path name, converts day of year 
to month and orders the file paths by time, where time = year + (month-0.5)/12. 

2. .hdf files are read and the NDVI variables for each .hdf file are extracted for a coodinate grid. 

3. Calculates mean NDVI per month over the time period

4. Calculates the anomaly for each NDVI reading by subtracting the monthly mean


INPUT;
MODIS HDF files containing the NDVI readings for each pixel location on a 7200x3600 coordinate grid

OUTPUT;
Lists for time series plots:
1.Monthly retrievals
2.Anomalies 
3.Precipitation readings 
4.Average monthly precipitation
5.Times when readings were taken


AUTHOR;
Freya Espir

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


def get_ndvi_data(file_outpath,path,suffix):
    F = (extract_files(path,'.hdf'))
    
    try:    
        file = open(file_outpath,'r')
    except IOError: # writing to file if it doesn't exist
        fileinfo = fetch_NDVI_file_attributes(F)
        file = open(file_outpath, 'w')
        with file as f:
            json.dump(fileinfo, f) 
            
    with open(file_outpath,'r') as reading:
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
    return YEAR,DAY,TIME,MONTH,FILE
       
def ndvi_data_for_coordinate_range(outpath,files,months,times,years,dataname,latbnds,lonbnds,scale_factor,add_offset):
    
    try:
        file = open(outpath,"r")
    except IOError:
        ndvi_means = []
        for item in files:
            readdata = get_NDVI_readings(item,latbnds,lonbnds,dataname,scale_factor,add_offset)
            mean_ndvi = grid_average(readdata)
            ndvi_means.append(mean_ndvi)
        with open(outpath, "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            for val in ndvi_means:
                writer.writerow([val])
           
    
    ndvi_vals = []
    with open(outpath,"r") as r:
        reader = csv.reader(r,delimiter=' ')
        for row in reader:   
            values = float((",".join(row)))
            ndvi_vals.append(values)
            
    NDVIs = np.asarray(ndvi_vals)
    
    seasonal_anoms = seasonal_anomalies(years,months,times,NDVIs,'NDVI yearly means time series dry season','NDVI yearly means time series wet season', 'NDVI means')
    #calculating monthly average and writing to list
    data = monthly_average(months,NDVIs)
    monthly_retrievals = data[1]
    NDVI_average = data[0]

    anomalies = [] 
    #calculating anomalies for each of the precipitation values
    for i in range(len(months)):
        anom = NDVIs[i] - NDVI_average[int(months[i])-1]
        anomalies.append(anom)

    #converting plotting list to array
    anomalies = np.asarray(anomalies)
    monthly_retrievals = np.asarray(monthly_retrievals)
    Anoms = np.array(anomalies)

    out = {'Monthly Retrievals':monthly_retrievals,'Anomalies': Anoms, 'NDVIs':NDVIs, 'Average NDVI':NDVI_average, 'Times':times, 'dry times': seasonal_anoms['dry times'], 'rainy times': seasonal_anoms['rainy times'], 'rainy anomalies': seasonal_anoms['rainy anomalies'], 'dry anomalies': seasonal_anoms['dry anomalies'],'nine year anoms wet': seasonal_anoms['nine year anomalies rain'],'nine year anoms dry': seasonal_anoms['nine year anomalies dry']}
    return out
    

def process_ndvi(lonbnds,latbnds,outpath,ndvi_path,suffix, ndvi_files_suffix, ndvi_means_suffix):
    print "Processing NDVI"
    ndvi_outpath = outpath + ndvi_files_suffix
    mean_ndvi_out = outpath + ndvi_means_suffix
    getdata = get_ndvi_data(ndvi_outpath,ndvi_path,suffix)
    attributes = ndvi_data_for_coordinate_range(mean_ndvi_out,getdata[4],getdata[3], getdata[2],getdata[0],'CMG 0.05 Deg Monthly NDVI',latbnds,lonbnds,10000.,0.)
    return attributes

 
def process_ndvi_mg(lonbnds,latbnds,outpath,ndvi_path,suffix):
    print "Processing NDVI"
    ndvi_outpath = outpath + "NDVI_files_mg.json"
    mean_ndvi_out = outpath + "mean_ndvi_out_mg.csv"
    getdata = get_ndvi_data(ndvi_outpath,ndvi_path,suffix)
    attributes = ndvi_data_for_coordinate_range(mean_ndvi_out,getdata[4],getdata[3], getdata[2],getdata[0],'CMG 0.05 Deg Monthly NDVI',latbnds,lonbnds,10000.,0.)
    return attributes