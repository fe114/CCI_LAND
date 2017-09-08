#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 09:26:49 2017

@author: fespir
"""

#cloud fraction 

import json
from numpy import * 
from math import *
import sys
import netCDF4 
from netCDF4 import Dataset
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
from operator import itemgetter 
from scipy import stats
from itertools import repeat
from scipy.stats import t 
from atsr import *
from file_search import *
from geolocation import *
import csv
from PAR_processing import *

def get_cloud_frac_files(outpath,path,suffix):
    try:    
        file = open(outpath,'r')
    except IOError: # writing to file if it doesn't exist
        cloud_frac = get_PAR_attributes(extract_files(path,suffix))
        file = open(outpath, 'w')
        with file as f:
            json.dump(cloud_frac, f) 

    with open(outpath, 'r') as f:
        attributes = json.load(f)
        dates  = attributes['Time']
        months = attributes['Month']
        years  = attributes['Year']
        files  = attributes['File']

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
    return TIMES,FILES,YEARS,MONTHS


#toa_swdn_allsky(lat,lon)
def read_cloud_frac_files(outpath,times,files,years,months,latname,lonname,variable,latval,lonval):
    cloud_frac_means = []   
 
    try: 
        file = open(outpath,"r")
    except IOError:
        for item in files:
            gettingdata = getdata_coordinategrid_lat_lon(item,lonname,latname,lonval,latval,variable)
            mean = grid_average(gettingdata)
            cloud_frac_means.append(float(mean))# getdata(filepath, longitude name in .nc file, latitude name in .nc file, variable name in .nc file, lon coordinate, lat coordinate)
        with open(outpath, "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            for val in cloud_frac_means:
                writer.writerow([val])  
            

                
    cloud_fracs = []
    with open(outpath,"r") as r:
        reader = csv.reader(r,delimiter=' ')
        for row in reader:   
            cloud_fracvals = float((",".join(row)))
            cloud_fracs.append(cloud_fracvals) 

    cloud_fracs = np.asarray(cloud_fracs)
    average = monthly_average(months,cloud_fracs)

    cloud_frac_monthly_average = average[0]
    Monthly_readings = average[1]

    
    #computing anomalie for each AOD value by substracting the monthly average
    anomalies = [] 

    for i in range(len(months)):
        anom = cloud_fracs[i] - cloud_frac_monthly_average[int(months[i])-1]
        anomalies.append(anom)

    seasonal_anoms = seasonal_anomalies(years,months,times,cloud_fracs,'Cloud Fraction yearly means time series dry season','Cloud Fraction yearly means time series wet season', 'Cloud Fraction means')
    
    
    #converting plotting list to array
    anomalies = np.asarray(anomalies)
    monthly_retrievals = np.asarray(Monthly_readings)
    out = {'Anomalies':anomalies, 'Monthly Retrievals':monthly_retrievals,'cloud_frac Monthly Average':cloud_frac_monthly_average, 'Times':times, 'cloud_fracs':cloud_fracs, 'dry times': seasonal_anoms['dry times'], 'rainy times': seasonal_anoms['rainy times'], 'rainy anomalies': seasonal_anoms['rainy anomalies'], 'dry anomalies': seasonal_anoms['dry anomalies'],'nine year anoms wet': seasonal_anoms['nine year anomalies rain'],'nine year anoms dry': seasonal_anoms['nine year anomalies dry']}
    return out

def process_cloud_frac(lat,lon,outpath,cloud_frac_path,suffix,Cloud_frac_fileinfo_suffix,Cloud_frac_data_suffix):
    cloud_frac_info_out = outpath + Cloud_frac_fileinfo_suffix
    cloud_frac_outpath = outpath + Cloud_frac_data_suffix
    getfile = get_cloud_frac_files(cloud_frac_outpath,cloud_frac_path,suffix)
    data = read_cloud_frac_files(cloud_frac_info_out,getfile[0],getfile[1],getfile[2],getfile[3],'lat', 'lon', 'cldfr_allsky',lat,lon)
    return data
"""
def process_cloud_frac_mg(outpath,cloud_frac_path,lat,lon,suffix):
    cloud_frac_info_out = outpath + "CCI_vals_mg.csv"
    CCI_outpath = outpath + "CCI_attributes_mg.json"
    getfile = get_CCI_files(CCI_outpath,CCI_path,suffix)
    data = read_CCI_files(CCI_info_out,getfile[0],getfile[1],getfile[2],getfile[3], 'lat', 'lon', 'cldfr_allsky',lat,lon)
    return data   
"""