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

def getdata(read,lon,lat,data,coord_lon,coord_lat):
    readfile = Dataset(read,mode='r')
    lons = readfile.variables[lon][:] 
    lats = readfile.variables[lat][:]
    datavalues = readfile.variables[data][:,:]
    print datavalues
    coordinate_lon = np.where(lons == coord_lon)
    coordinate_lat = np.where(lats == coord_lat)
    data_at_point = datavalues[coordinate_lat,coordinate_lon]
    
   #data_at_point = float(data_at_point)
    #dataout =  {"AATSR_READING" : data_at_point, "LON_LOCATION" : coordinate_lon, "LAT_LOCAION": coordinate_lat}
    dataout = data_at_point
    return dataout


def array2json(your_array):
    your_array = pd.Series(your_array).to_json(orient='values')
    return your_array

def remove_nan(list_unfiltered,times,months,years,filename):
    nonan = ~np.isnan(list_unfiltered)
    y = list_unfiltered[nonan]
    x = times[nonan]
    MONTHS = months[nonan]
    YEARS = years[nonan]
    FILES = filename[nonan]    
    outlists = {"AOD": y, "TIMES": x,"MONTHS":MONTHS, "YEARS" : YEARS, "FILES" : FILES}
    
    return outlists