
"""
Created on Mon Aug 14 14:49:38 2017

@author: fespir
# Name: General moduel used by all data sets for processing lat lon location

# Functions: get_data, array2json,  remove_nan,get_coordinate_grid_lat_lon, get_coordinate_grid_lon_lat, grid_average, monthly_average



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
import json
import pandas as pd
import pickle

#reads a .nc file and uses the lat and lon location for a single point and outputs the data reading for that point
def getdata(read,lat,lon,data,coord_lat,coord_lon):
    readfile = Dataset(read,mode='r')
    lons = readfile.variables[lon][:] 
    lats = readfile.variables[lat][:]
    datavalues = readfile.variables[data][:,:]
    coordinate_lon = np.where(lons == coord_lon)
    coordinate_lat = np.where(lats == coord_lat)
    data_at_point = datavalues[coordinate_lat,coordinate_lon]
    data_at_point = float(data_at_point)
    return data_at_point
#converts an array to json format for savefile
def array2json(your_array):
    your_array = pd.Series(your_array).to_json(orient='values')
    return your_array

#removes nan values from a list 
#sorts through the equivalent y,x,months,files,years lists and removed the data for the nan index vales
def remove_nan(list_unfiltered,times,months,years,filename):
    nonan = ~np.isnan(list_unfiltered)
    y =list_unfiltered[nonan]
    x = times[nonan]
    MONTHS = months[nonan]
    YEARS = years[nonan]
    FILES = filename[nonan]
    return y,x,MONTHS,YEARS,FILES

#reads .nc file and extracts data for a coordinate GRID. note: input as lon:lat                           
def getdata_coordinategrid_lon_lat(read,lon,lat,lon_boundary,lat_boundary,variable):
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
#reads .nc file and extracts data for a coordinate GRID. note: input as lat:lon  
def getdata_coordinategrid_lat_lon(read,lon,lat,lon_boundary,lat_boundary,variable):
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
        precip = readfile.variables[variable][lat_zero:lat_max,lon_zero:lon_max]
        return precip   

#calculates the mean value for a coodinate grid  
# input = data from .nc file  
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
  
#calculates the monthly average  
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