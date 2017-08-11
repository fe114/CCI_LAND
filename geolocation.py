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
import json
import pandas as pd
import pickle

def getdata(read,lon,lat,data,coord_lon,coord_lat):
    readfile = Dataset(read,mode='r')
    lons = readfile.variables[lon][:] 
    lats = readfile.variables[lat][:]
    datavalues = readfile.variables[data][:,:]
    coordinate_lon = np.where(lons == coord_lon)
    coordinate_lat = np.where(lats == coord_lat)
    data_at_point = datavalues[coordinate_lat,coordinate_lon]
    data_at_point = float(data_at_point)
    return data_at_point


def array2json(your_array):
    your_array = pd.Series(your_array).to_json(orient='values')
    return your_array

def remove_nan(list_unfiltered,times,months,years,filename):
    nonan = ~np.isnan(list_unfiltered)
    print 'no_nan', nonan
    y =list_unfiltered[nonan]
    x = times[nonan]
    MONTHS = months[nonan]
    YEARS = years[nonan]
    FILES = filename[nonan]
    return y,x,MONTHS,YEARS,FILES
                            
def getdata_coordinategrid(read,lon,lat,lat_boundary,lon_boundary,variable):
        readfile = Dataset(read,mode='r')
        #print readfile.variables.keys()
        lons = readfile.variables[lon][:]
        lats = readfile.variables[lat][:]
        lat_inds = np.where((lats > lat_boundary[0]) & (lats < lat_boundary[1]))[0]
        lon_inds = np.where((lons > lon_boundary[0]) & (lons < lon_boundary[1]))[0]
        lat_zero = lat_inds[0]
        lat_max = lat_inds[len(lat_inds)-1]
        lat_zero, lat_max
        lon_zero = lon_inds[0]
        lon_max = lon_inds[len(lon_inds)-1]
        precip = readfile.variables[variable][lon_zero:lon_max,lat_zero:lat_max]
        return precip
    
def grid_average(precip):
        sum_list = []
        ncols = len(precip[0])
        for col in range(ncols): #add each row in NDVI
            sum_list.append(sum(row[col] for row in precip))
        #get the total mean ndvi for lat lon range 
        sum_precip = sum(sum_list) #add the summed rows to get one total mean ndvi value
        N_values = precip.size # number of NDVI readings in lat,lon range
        mean_precip = sum_precip/N_values # mean NDVI value for the coordinate range
        #all_precip_means.append(mean_precip)
        return mean_precip
    
def monthly_average(months,data):
    count = 1
    month_list = linspace(1,12,12)
    monthly_averages = []
    monthly_readings = []#calculating monthly average
    for i in range(len(month_list)):
        location = np.where(months == count)
        location = np.asarray(location)
        monthly_precipitation = data[location]
        total = sum(monthly_precipitation)
        length = monthly_precipitation.size
        average = total/length
        monthly_averages.append(average)
        monthly_readings.append(length)
        count = count+1 
    return monthly_averages, monthly_readings