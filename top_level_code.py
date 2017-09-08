#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
NAME;
Top level programme to process data using sub processing scripts for a cordinate grid using lat and lon boundaries



PURPOSE;
1. Program uses sub processing scripts to process data for one lat lon cooridnate range determined from the raw_input function 
 
2. Dictionaries are produced for each of the following datasets:
AOD, Precipitation, NDVI, PAR_diffuse,PAR_total, toa_allsky, toa_par_tot_allsky and Cloud Fraction 

3. Dictionaries are saved to a .txt file using the pickle function. The name will be catagorised by the
lat and lon location followed by the data name, e.g. for Aerosol the file name is '-10-64precip.txt'
for the lat lon boundary -11,-9 N -65,-63E.

4. Use the data2plot.py program to plot data from the files

AUTHOR;
Freya Espir

'''

from numpy import * 
import math
from math import *
from netCDF4 import Dataset
import os,sys
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt 
from scipy import stats 
from atsr import *
from file_search import *
from geolocation import *
from AOD import *
from precip import *
from ndvi import *
from PAR import * 
from Cloud_Frac import *
from pltplots import *
import csv
import pickle
from matplotlib import gridspec
from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go
from scipy.stats import linregress

figpath = '/group_workspaces/cems2/nceo_generic/CCI_LAND/figs/rondonia/'
aerosol_path = "/group_workspaces/cems/aerosol_cci/public/cci_products/AATSR_ORAC_v04-01/L3_MONTHLY/"
precip_path = "/group_workspaces/cems2/nceo_generic/satellite_data/precipitation/TRMM/monthly"
ndvi_path = "/group_workspaces/cems2/nceo_generic/satellite_data/modis_c6/ndvi/"
outpath = "/group_workspaces/cems2/nceo_generic/CCI_LAND/output_data/"
PAR_path = "/group_workspaces/cems/cloud_ecv/public/ESA_Cloud_CCI/CLD_PRODUCTS/FLX/L3C/AATSR_ENVISAT/v2.1/"
cloud_frac_path = "/group_workspaces/cems/cloud_ecv/public/ESA_Cloud_CCI/CLD_PRODUCTS/FLX/L3C/AATSR_ENVISAT/v2.0/"


#lat and lon bounds for region which hasn't undergone deforestation
#lat_bnds_controlregion, lon_bnds_controlregion = [-9,-7], [-67,-65]

#lat and lon bounds rondonia
#lat_bnds, lon_bnds = [-11,-9], [-65,-63]

lat_bnds = []
lon_bnds = []

lat1 = lat_bnds.append(int(raw_input('lat lower boundary ')))
lat2 = lat_bnds.append(int(raw_input('lat upper boundary ')))
lon1 = lon_bnds.append(int(raw_input('lon lower boundary ')))
lon2 = lon_bnds.append(int(raw_input('lon upper boundary ')))

lat = str(sum(lat_bnds)/len(lat_bnds))
lon = str(sum(lon_bnds)/len(lon_bnds))

#processing data for Rondonia where rapid deforestation has occured since ~ 1990 

def dicttofile(filename,data):
    with open(filename, 'wb') as handle:
        pickle.dump(data, handle)
        
Aerosol_data = process_aerosol(lon_bnds,lat_bnds,outpath,aerosol_path,lat + lon +'aodfiles.json',lat + lon +'aodlonlat.csv',lat + lon +'AOD_no_nan.json')
dicttofile(outpath+lat+lon+'aod.txt',Aerosol_data)

Precipitation_data = process_precip(lon_bnds,lat_bnds,outpath,precip_path,lat + lon +'precipitation.txt', lat + lon +'precip_means.csv') 
dicttofile(outpath+lat+lon+'precip.txt', Precipitation_data)

NDVI_data = process_ndvi(lon_bnds,lat_bnds,outpath,ndvi_path,".hdf",lat + lon +"NDVI_files.json", lat + lon +"mean_ndvi_out.csv")
dicttofile(outpath+lat+lon+'ndvi.txt', NDVI_data)

PAR_Total_data = process_PAR(outpath,PAR_path,lat_bnds,lon_bnds,"boa_par_tot_allsky","fv2.1.nc",lat + lon +"PAR_total_vals.csv",lat + lon +"par_total_attributes.json") #base of atmosphere total allsky composite w/m^2
dicttofile(outpath+lat+lon+'PAR_total.txt',PAR_Total_data)

PAR_Diffuse_data = process_PAR(outpath,PAR_path,lat_bnds,lon_bnds,"boa_par_dif_allsky", "fv2.1.nc", lat + lon +"PAR_diffuse_vals.csv",lat + lon +"par_diffuse_attributes.json" )#base of atmosphere diffuse allsky composite w/m^2
dicttofile(outpath+lat+lon+'PAR_diffuse.txt', PAR_Diffuse_data)

cloud_frac_data = process_cloud_frac(lat_bnds,lon_bnds,outpath, cloud_frac_path, "v2.0.nc",lat + lon +"cloud_frac_vals.csv", lat + lon +"cloud_frac_attributes.json")
dicttofile(outpath+lat+lon+ 'cloudfrac.txt', cloud_frac_data)

toa_allsky_data = process_PAR(outpath,PAR_path,lat_bnds, lon_bnds, "toa_swdn_allsky","fv2.1.nc", lat + lon +"toa_allsky_vals.csv",lat + lon + "toa_attributes.json")
dicttofile(outpath+lat+lon+ 'toa_allsky.txt', toa_allsky_data)

toa_par_tot_allsky_data = process_PAR(outpath,PAR_path,lat_bnds, lon_bnds, "toa_par_tot_allsky","fv2.1.nc", lat + lon +"toa_par_tot_allsky_vals.csv",lat + lon + "toa_par_tot_allsky_attributes.json" )
dicttofile(outpath+lat+lon+'toa_par_tot_allsky.txt',toa_par_tot_allsky_data)

#--------------------------------------------
