#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
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
Lists for time series plots:
Monthly retrievals
Anomalies 
Precipitation readings 
Average monthly precipitation
Times when readings were taken


AUTHOR;
Freya Espir

"""


import json
from numpy import * 
from math import *
from netCDF4 import Dataset
import numpy as np
from atsr import *
from file_search import *
from geolocation import *
from AOD import *
import csv

def get_precip_data(file_outpath,path,suffix):
    try:    
        file = open(file_outpath,'r')
    except IOError: # writing to file if it doesn't exist
        precipitation = fetch_aerosol_file_attributes(extract_files(path,suffix))
        file = open(file_outpath, 'w')
        with file as f:
            json.dump(precipitation, f) 

    # Reading data back
    with open(file_outpath, 'r') as f:
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
    MONTHS = m[SORT] 
    YEARS = yr[SORT]
    FILES = f[SORT]
    
    return SORT,TIMES,MONTHS,YEARS,FILES
                         
def data_for_coordinate_range(outpath,files,times,latname,lonname,dataname,latbnds,lonbnds,months):
    try: 
        file = open(outpath, 'r')
    except IOError:
        all_precip_means = []
        for item in files:
            precip = getdata_coordinategrid(item, lonname,latname,latbnds,lonbnds,dataname)                       
            mean_precip = grid_average(precip)
            all_precip_means.append(mean_precip)
        with open(outpath, "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            for val in all_precip_means:
                writer.writerow([val]) 

#ndvi = outpath + 'ndvi_mean.csv'

    precip_vals = []
    with open(outpath,"r") as r:
        reader = csv.reader(r,delimiter=' ')
        for row in reader:   
            values = float((",".join(row)))
            precip_vals.append(values)
    precipitations = np.asarray(precip_vals)


    #calculating monthly average and writing to list
    data = monthly_average(months,precipitations)
    monthly_retrievals = data[1]
    precip_average = data[0]

    anomalies = [] 

    #calculating anomalies for each of the precipitation values
    for i in range(len(months)):
        anom = precipitations[i] - precip_average[int(months[i])-1]
        anomalies.append(anom)

    #converting plotting list to array
    anomalies = np.asarray(anomalies)
    monthly_retrievals = np.asarray(monthly_retrievals)
    Anoms = np.array(anomalies)
    out = {'Monthly Retrievals':monthly_retrievals,'Anomalies': Anoms, 'Precipitation':precipitations, 'Average Precipitation':precip_average, 'Times':times}
    return out
    
def process_precip(lonbnds,latbnds,outpath,precip_path):
    precip_outpath = outpath+'precipitation.txt'
    mean_precip_out = outpath = 'precip_means.csv'
    getdata=get_precip_data(precip_outpath,precip_path,"_TRMM_3B43.7.nc")
    data_coord_range = data_for_coordinate_range(mean_precip_out,getdata[4],getdata[1],'nlat','nlon','precipitation',latbnds,lonbnds,getdata[2])
    
    return data_coord_range