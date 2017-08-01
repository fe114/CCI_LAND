#!/usr/bin/env python2
# -*- coding: utf-8 -*-


"""
NAME;
Plot of the global distribution of the mean cloud fraction for June 2008.

PURPOSE;
This program to processes CFC, lon and lat variables from .nc file and outputs a figure plot of
the global distribution of mean cloud fraction.

DESCRIPTION;
The processor only requires one AATSR dataset

A series of functions are applied to retrieve variable information from 
AATSR dataset and change the projection geometry. 

The output figure is saved as a .pdf file

INPUT;
NetCDF file containing the CFC variables for each pixel location

OUTPUT;
PDF file of the figure plot.

AUTHOR;
Freya Espir
"""

from numpy import * 
from math import *
import os,sys
from netCDF4 import Dataset
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#paths
figpath = '/group_workspaces/cems2/nceo_generic/CCI_LAND/figs/'
filein = '/group_workspaces/cems/cloud_ecv/public/ESA_Cloud_CCI/CLD_PRODUCTS/L3C/AATSR_ENVISAT/v2.0/2008/200806-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AATSR_ENVISAT-fv2.0.nc'

#read in lon, lat, cfc data from .nc file
ncfile = Dataset(filein, mode='r') 
lons = ncfile.variables['lon'][:] #retreive longitude variables 
lats = ncfile.variables['lat'][:] #atitude variables
cfc = ncfile.variables['cfc'][0,:,:] #cfc variables. (time,lat,lon)
ncfile.close() 
#print the file format;
print ncfile.file_format 

#change size of figure plot(lon,lat)
fig = plt.figure(figsize=(40,30)) 

#Set the GeoAxes to a 'standard' projection onto a cylinder tangent at the Equator
ax = plt.axes(projection=ccrs.PlateCarree())

#change_geometry(numrows, numcols, numplots)
ax.change_geometry(3,2,1)

#labeling the gridlines with numbers, colour line black with the style -- and a linewidth of 2. 
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                 linewidth=2, color='black', alpha=0.5, linestyle='--')

#marking on the longitude gridlines.
x = linspace(-180, 180, 13) # (start,stop,number of increments)
#mark longitude on x axis
gl.xlocator = mticker.FixedLocator(x) 

#same for latitude 
y = linspace(-90,90,13)
gl.ylocator = mticker.FixedLocator(y)

#removing lon and lat from top and left axis
gl.xlabels_top = False
gl.ylabels_left = False

#create contour plot using the lat,lon and cfc variables.
plt.contourf(lons, lats, cfc, 5,
             transform=ccrs.PlateCarree())

#Add colourbar. Shrink indicates the fraction by which to shrink the colourbar
plt.colorbar(ax=ax, shrink=.8)

#x.formatter makes the -ve x values west and positive values east
#y.formatter -ve = south, + = north
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

#mark on coastlines 
ax.coastlines()
plt.title('Cloud Fraction June 2008')
#show the plot
plt.show()

plt.savefig(figpath + 'cfc_june_08.pdf')
#This is a good webpage giving information to do different plots
#http://wrf-python.readthedocs.io/en/latest/plot.html