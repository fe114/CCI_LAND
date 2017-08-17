#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 14:49:38 2017

@author: fespir
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

def single_location(f_out,FILES,lon,lat):
    try: 
        file = open(f_out, 'r')
    except IOError:
        aodmeans = []
        for item in FILES:
            gettingdata = getdata(item,'longitude', 'latitude', 'AOD550_mean',lon,lat)
            aodmeans.append(gettingdata)# getdata(filepath, longitude name in .nc file, latitude name in .nc file, variable name in .nc file, lon coordinate, lat coordinate)
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

def removing_nans(no_nans,AOD_VALUE,SORT,TIMES,MONTHS,YEARS,FILES,outpath):  
    #removing item from list if the corresponding AOD value is = nan
    
    A = AOD_VALUE[no_nans]
    S = SORT[no_nans]
    T = TIMES[no_nans]
    M = MONTHS[no_nans]
    Y = YEARS[no_nans]
    F = FILES[no_nans]
    
    #converting A array to json format
    Y = array2json(A)
 
    try: 
        file = open(outpath, 'r')
    except IOError:
        file = open(outpath, 'w')
        with file as towrite:
            json.dump(Y, towrite)
    
    with open(outpath, 'r') as toread:
        data = json.load(toread)

    
    month_list = linspace(1,12,12)

    #Calculating monthly AOD average
    count = 1
    AOD_average = []
    monthly_retrievals = []
    
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

    for i in range(len(M)):
        anom = A[i] - AOD_average[int(M[i])-1]
        anomalies.append(anom)


    #converting plotting list to array
    anomalies = np.asarray(anomalies)
    monthly_retrievals = np.asarray(monthly_retrievals)
    out = {'Anomalies':anomalies , 'Monthly Retrievals':monthly_retrievals,'AOD Monthly Average':AOD_average, 'Times':T, 'AOD Measurements':A}
    return out

def process_aerosol(lon,lat,outpath,aerosol_path):
    aerfile = outpath+'aodfiles.json'
    AOD_out = outpath + 'aodlonlat.csv'
    AOD_no_nan= outpath + 'AOD_no_nan.json'
    Aerosol_data = get_aerosol_data(aerfile,aerosol_path,"04.01.nc")
    Aerosol_one_coordinate = single_location(AOD_out,Aerosol_data[4],lon,lat)
    Remove_nan = removing_nans(Aerosol_one_coordinate[0],Aerosol_one_coordinate[1],Aerosol_data[0],Aerosol_data[1],Aerosol_data[2],Aerosol_data[3],Aerosol_data[4],AOD_no_nan)
    return Remove_nan
   
