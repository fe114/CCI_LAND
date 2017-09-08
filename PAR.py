#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 13:14:04 2017

@author: fespir
"""
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
from solar import *

def get_PAR_files(outpath,path,suffix):
   
    try:    
        file = open(outpath,'r')
    except IOError: # writing to file if it doesn't exist
        PAR = get_PAR_attributes(extract_files(path,suffix))
        file = open(outpath, 'w')
        with file as f:
            json.dump(PAR, f) 
    
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
def read_PAR_files(outpath,files,months,times,years,latname,lonname,variable,latval,lonval):
    PAR_means = []   
    
    try: 
        file = open(outpath,"r")
    except IOError:
        for item in files:
            gettingdata = getdata_coordinategrid_lat_lon(item,lonname,latname,lonval,latval,variable)
            mean = grid_average(gettingdata)
            PAR_means.append(float(mean))# getdata(filepath, longitude name in .nc file, latitude name in .nc file, variable name in .nc file, lon coordinate, lat coordinate)
        with open(outpath, "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            for val in PAR_means:
                writer.writerow([val])  
         
    PARs = []
    with open(outpath,"r") as r:
        reader = csv.reader(r,delimiter=' ')
        for row in reader:   
            PARvals = float((",".join(row)))
            PARs.append(PARvals) 

    PARs = np.asarray(PARs)
    latitude = ((latval[0] + latval[1])/2.0)
    
    solar_rad_daily = []
    for i in range(len(times)):
        caldat = ((times[i] %1.0)* 365)
        solar_rad_daily.append(solar_radiation_daily(latitude,caldat))
    
    
    
    #print PAR_correction, PARs
    
    seasonal_anoms = seasonal_anomalies(years,months,times,PARs,'PAR yearly means time series dry season','PAR yearly means time series wet season', 'PAR means')
    
    average = monthly_average(months,PARs)

    PAR_monthly_average = average[0]
    Monthly_readings = average[1]

    
    #computing anomalie for each AOD value by substracting the monthly average
    anomalies = [] 
    standardised_anoms = []
    st_d = np.std(PAR_monthly_average)
    
    
    for i in range(len(months)):
        anom = PARs[i] - PAR_monthly_average[int(months[i])-1]
        standardised_anoms.append(anom/st_d)
        anomalies.append(anom)

    
    #converting plotting list to array
    anomalies = np.asarray(anomalies)
    monthly_retrievals = np.asarray(Monthly_readings)
    out = {'Anomalies':anomalies, 'Monthly Retrievals':monthly_retrievals,'PAR Monthly Average':PAR_monthly_average, 'Times':times, 'PARs':PARs, 'dry times': seasonal_anoms['dry times'], 'rainy times': seasonal_anoms['rainy times'], 'rainy anomalies': seasonal_anoms['rainy anomalies'], 'dry anomalies': seasonal_anoms['dry anomalies'], 'PAR monthly data': PARs, 'Standardised Anomalies': standardised_anoms, 'Daily Solar Radiation' : solar_rad_daily}
    return out

def process_PAR(outpath,PAR_path,lat,lon,variable,suffix,PAR_file_info_suffix, PAR_data_suffix):
    PAR_info_out = outpath + PAR_file_info_suffix
    PAR_outpath = outpath + PAR_data_suffix
    getfile = get_PAR_files(PAR_outpath,PAR_path,suffix)
    data = read_PAR_files(PAR_info_out,getfile[1],getfile[3],getfile[0],getfile[2],'lat', 'lon', variable,lat,lon)
    return data
  
"""
def process_PAR_total(outpath,PAR_path,lat,lon,variable,suffix):
    PAR_info_out = outpath + "PAR_total_vals.csv"
    PAR_outpath = outpath + "par_total_attributes.json"
    getfile = get_PAR_files(PAR_outpath,PAR_path,suffix)
    data = read_PAR_files(PAR_info_out,getfile[1],getfile[3],getfile[0],getfile[2],'lat', 'lon', variable,lat,lon)
    return data

def process_PAR_total_mg(outpath,PAR_path,lat,lon,variable,suffix):
    PAR_info_out = outpath + "PAR_total_vals.csv"
    PAR_outpath = outpath + "par_total_attributes.json"
    getfile = get_PAR_files(PAR_outpath,PAR_path,suffix)
    data = read_PAR_files(PAR_info_out,getfile[1],getfile[3],getfile[0],getfile[2],'lat', 'lon', variable,lat,lon)
    return data

def process_PAR_diffuse_mg(outpath,PAR_path,lat,lon,variable,suffix):
    PAR_info_out = outpath + "PAR_diffuse_mg_vals.csv"
    PAR_outpath = outpath + "par_diffuse_mg_attributes.json"
    getfile = get_PAR_files(PAR_outpath,PAR_path,suffix)
    data = read_PAR_files(PAR_info_out,getfile[1],getfile[3],getfile[0],getfile[2],'lat', 'lon', variable,lat,lon)
    return data
"""