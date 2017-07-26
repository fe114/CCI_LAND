#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 10:33:15 2017

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

print lats
#AOD_mean_units = readfile.variables['AOD550_mean'].units


#change the size of the figure plot (lon,lat)
fig = plt.figure(figsize=(30,40)) 


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

#removing labels from the top and left axes
gl.xlabels_top = False
gl.ylabels_left = False

#creating the contour plot using the lat,lon and AOD_mean variables.
levels = linspace(0,1,100)
plot = plt.contourf(lons, lats, AOD_mean, levels, 
             transform=ccrs.PlateCarree())

#create a colourbar. Shrink indicates the fraction by which to shrink the colourbar
#colourbar = plt.colorbar(plot)
plt.colorbar(ax=ax, shrink=.4)

#x.formatter makes the -ve x values west and positive values east
#y.formatter -ve = south, + = north
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

#mark on coastlines 
ax.coastlines()
plt.title('Global distribution of Aerosol Optical Depth at 550nm for June 2008')
plt.savefig( figpath + 'AOD550_june_2008.pdf')

#show the plot
plt.show()

coordlon = np.where(lons == -63.5)
coordlat = np.where(lats == -11.5)
print AOD_mean[coordlon,coordlat]
