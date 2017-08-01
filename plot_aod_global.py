#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 17:03:59 2017

@author: fespir
"""

"""
NAME;
Plot of the global distribution of Aerosol Optical Depth at 550nm for June 2008

PURPOSE;
This program to processes AOD, lon and lat variables for each pixel from a
.nc file taken by AATSR and produces a figure plot of the AOD data with matplotlib
using a colour spectrum with 100 levels. 

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
import xarray as xr
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.backends.backend_pdf import PdfPages

figpath = '/group_workspaces/cems2/nceo_generic/CCI_LAND/figs/'

#retrieve data from the folder. I saved the .nc file to the same folder as the one this .py script is saved to
#so I don't need to located it
Aero_data = '/group_workspaces/cems/aerosol_cci/public/cci_products/AATSR_ORAC_v04-01/L3_MONTHLY/2008/200806-ESACCI-L3C_AEROSOL-AER_PRODUCTS-AATSR-ENVISAT-ORAC-MONTHLY-fv04.01.nc'
readfile = Dataset(Aero_data, mode='r') #reading the data
print readfile.file_format
#print readfile.variables['AOD550_mean']
lons = readfile.variables['longitude'][:] #retreive the longitude variables from the data file
lats = readfile.variables['latitude'][:] #for the latitude variables
AOD_mean = readfile.variables['AOD550_mean'][:,:]  
readfile.close() 
print "AOD_mean =\n",AOD_mean[:,:]
print "latitudes =\n",lats[:]
print "longitudes =\n",lats[:]
#get the units for AOD_mean
#AOD_mean_units = readfile.variables['AOD550_mean'].units


#change the size of the figure plot (lon,lat)
fig = plt.figure(figsize=(40,30)) 
#fig = plt.figure(figsize=(80,60)) 

#Set the GeoAxes to a projection onto a cylinder tangent at the Equator.
ax = plt.axes(projection=ccrs.PlateCarree())

#change_geometry(numrows, numcols, num)
ax.change_geometry(3,2,1)

#labeling the gridlines with just numbers, colouring the line black with the style -- and a linewidth of 2. 
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                 linewidth=2, color='black', alpha=0.5, linestyle='--')

#marking on the longitude gridlines. Linspace(start,stop,number of increments)
x = linspace(-180, 180, 13)
#print x 
gl.xlocator = mticker.FixedLocator(x) #mark longitude on x axis

#same for latitude 
y = linspace(-90,90,13)
gl.ylocator = mticker.FixedLocator(y)


gl.xlabels_top = False
gl.ylabels_left = False

#creating the contour plot using the lat,lon and AOD_mean variables.
levels = linspace(0,1,100)
plot = plt.contourf(lons, lats, AOD_mean, levels, 
             transform=ccrs.PlateCarree())

#create a colourbar. Shrink indicates the fraction by which to shrink the colourbar

#plt.pcolor(lons, lats, AOD_mean, cmap='cool')
colourbar = plt.colorbar(plot)
#x.formatter makes the -ve x values west and positive values east
#y.formatter -ve = south, + = north
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

#mark on coastlines 
ax.coastlines()
plt.title('Global distribution of Aerosol Optical Depth at 550nm for June 2008')
plt.savefig( figpath + 'aerosol.pdf')

#show the plot
plt.show()

coordlon = np.where(lons == 74.12)
coordlat = np.where(lats == 15.5)
print coordlon,coordlat