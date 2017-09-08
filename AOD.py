#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 14:49:38 2017

@author: fespir
# Name: Module to process aerosol data
# Functions: get_aerosol_data, single_location, removing_nans, process_aerosol
# Purpose: script recieves input from top level code (lat boundary, lon boundary, outpath and filenames) and outputs a dictionary of AOD data
# Output: 
'Anomalies', 'Monthly Retrievals','AOD Monthly Average', 'Times','AOD Measurements', 'dry times', 
'rainy times','rainy anomalies', 'dry anomalies', 'nine year anoms wet','nine year anoms dry','Standardised Anomalies'


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
import sys
import csv


#finds aserosol files and sorts files by the filepath name into date order
def get_aerosol_data(file_outpath,path,suffix):

    try: 
        file = open(file_outpath, 'r')
    except IOError:  
        file = open(file_outpath, 'w')
        aerosol_files = fetch_aerosol_file_attributes(extract_files(path,suffix))
        with file as f:
            json.dump(aerosol_files, f)


    # Reading data back
    with open(file_outpath, 'r') as f:
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
    
    return SORT,TIMES,MONTHS,YEARS,FILES

#uses sorted filenames to extract aod data for lat lon region
def single_location(f_out,FILES,lon,lat):
    
    try: 
        file = open(f_out, 'r')
    except IOError:
        aodmeans = []
        for item in FILES: # read,lon,lat,lat_boundary,lon_boundary,variable
            gettingdata = getdata_coordinategrid_lat_lon(item,'longitude','latitude',lon,lat,'AOD550_mean')
            mean = grid_average(gettingdata)
            aodmeans.append(float(mean))# getdata(filepath, longitude name in .nc file, latitude name in .nc file, variable name in .nc file, lon coordinate, lat coordinate)
        with open(f_out, "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            for val in aodmeans:
                writer.writerow([val])
   
    AODs = []
    with open(f_out,"r") as r:
        reader = csv.reader(r,delimiter=' ')
        for row in reader: 
            AODvals = float((",".join(row)))
            AODs.append(AODvals)    
    
    AOD_VALUE = np.asarray(AODs)
    #removing nan values
    #~np.isnan finds index of all non nan values in AODs list
    nonan = np.where(~np.isnan(AODs))
    return nonan, AOD_VALUE
#removes the 'nan' values within lat lon region and wites output to file

def removing_nans(no_nans,AOD_VALUE,SORT,TIMES,MONTHS,YEARS,FILES,outpath):  
    #removing item from list if the corresponding AOD value is = nan
    
    A = AOD_VALUE[no_nans]
    S = SORT[no_nans]
    T = TIMES[no_nans]
    M = MONTHS[no_nans]
    Y = YEARS[no_nans]
    F = FILES[no_nans]
    

    
    seasonal_anoms = seasonal_anomalies(Y,M,T,A, 'AOD yearly means time series dry seasnon','AOD yearly means time series wet season', 'AOD means')
    
  
        
    #converting A array to json format
    Anoms = array2json(A)
    
    try: 
        file = open(outpath, 'r')
    except IOError:
        file = open(outpath, 'w')
        with file as towrite:
            json.dump(Anoms, towrite)
  
    
    month_list = linspace(1,12,12)

    #Calculating monthly AOD average
    count = 1
    AOD_average = []
    monthly_retrievals = []
    
    #caulculating one average AOD value for each month
   
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
    standardised_anoms = []
    st_d = np.std(AOD_average)
    
    for i in range(len(M)):
        anom = A[i] - AOD_average[int(M[i])-1]
        standardised_anoms.append(anom/st_d)
        anomalies.append(anom)
    
    
    
    #converting plotting list to array
    anomalies = np.asarray(anomalies)
    monthly_retrievals = np.asarray(monthly_retrievals)
    
    
    
    out = {'Anomalies':anomalies, 'Monthly Retrievals':monthly_retrievals,'AOD Monthly Average':AOD_average, 'Times':T, 'AOD Measurements':A, 'dry times': seasonal_anoms['dry times'], 'rainy times': seasonal_anoms['rainy times'], 'rainy anomalies': seasonal_anoms['rainy anomalies'], 'dry anomalies': seasonal_anoms['dry anomalies'], 'nine year anoms wet': seasonal_anoms['nine year anomalies rain'],'nine year anoms dry': seasonal_anoms['nine year anomalies dry'], 'Standardised Anomalies': standardised_anoms}
    return out

def process_aerosol(lon_bnds,lat_bnds,outpath,aerosol_path, aerfile_suffix,AOD_out_suffix, AOD_no_nan_suffix):
    print "Processing Aerosol"
    aerfile = outpath+ aerfile_suffix
    AOD_out = outpath + AOD_out_suffix
    AOD_no_nan= outpath + AOD_no_nan_suffix
    Aerosol_data = get_aerosol_data(aerfile,aerosol_path,"04.01.nc")
    Aerosol_one_coordinate = single_location(AOD_out,Aerosol_data[4],lon_bnds,lat_bnds)
    Remove_nan = removing_nans(Aerosol_one_coordinate[0],Aerosol_one_coordinate[1],Aerosol_data[0],Aerosol_data[1],Aerosol_data[2],Aerosol_data[3],Aerosol_data[4],AOD_no_nan)
    return Remove_nan

