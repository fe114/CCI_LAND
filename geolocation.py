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

def getdata(point):
    readfile = Dataset(point,mode='r') #read data 
    lons = readfile.variables['longitude'][:] 
    lats = readfile.variables['latitude'][:]
    mean_AOD_values = readfile.variables['AOD550_mean'][:,:]
    coordlon = np.where(lons == -63.5)
    coordlat = np.where(lats == -11.5)
    AOD_at_point = mean_AOD_values[coordlat,coordlon]
    aod_out =  {"AOD_VAL" : AOD_at_point, "LON_LOCATION" : coordlon, "LAT_LOCAION": coordlat} 
    return aod_out