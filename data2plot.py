#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
#---------------------------------------------------------------------------------------------
# Name: Module to read and plot data
# Functions: datafile, tolog
# Purpose: script reads data files and produces various plots. Some include best fit lines & tstatistics 

1. After processing data in top_level_code.py, input lat and lon coodrdinate boundarys when prompted by the raw_input function
2. Script reads files containing one python dictionary for each data type
3. data types include AOD, Cloud Fraction, PAR_diffuse, PAR_total, PAR_direct, toa_allsky, NDVI, Precipitation, toa_par_tot_allsky
4. extracts lists from the dictionaries for plotting
plot types include:
    scatter, 3x1 plots with shared x axis, 3x2 plots with shared x axis, one axis plot with 3 datasets 
5. saves figures as a .pdf to 'figpath' location. 
6. The figure name will be catagorised by the lat an lon location 

#---------------------------------------------------------------------------------------------


"""
Created on Thu Sep  7 17:00:34 2017

@author: fespir
"""
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
outpath = "/group_workspaces/cems2/nceo_generic/CCI_LAND/output_data/"
'''
---------------------
function to read text files containing a python dictionary
for AOD, the datafile will a dictionary containing the return the following lists:
'Anomalies','Monthly Retrievals','AOD Monthly Average','Times','AOD Measurements','dry times', 'rainy times', 'rainy anomalies', 'dry anomalies', 'nine year anoms wet','nine year anoms dry', 'Standardised Anomalies'
for example to extract AOD monthly readings from the dictionary:
    AOD = datafile(outpath,lat,lon,AODsuffix)['AOD Measurements']
--------------------
'''

def datafile(outpath,lat,lon,suffix):
    with open(outpath+lat+lon+suffix, 'rb') as handle:
        dictionary = pickle.loads(handle.read())
        return dictionary

lat_bnds = []
lon_bnds = []

lat1 = lat_bnds.append(int(raw_input('lat lower boundary ')))
lat2 = lat_bnds.append(int(raw_input('lat upper boundary ')))
lon1 = lon_bnds.append(int(raw_input('lon lower boundary ')))
lon2 = lon_bnds.append(int(raw_input('lon upper boundary ')))

lat = str(sum(lat_bnds)/len(lat_bnds))
lon = str(sum(lon_bnds)/len(lon_bnds))

#reading the files
AOD_dict = datafile(outpath,lat,lon,'aod.txt')
PRECIP_dict = datafile(outpath,lat,lon,'precip.txt')
NDVI_dict = datafile(outpath,lat,lon,'ndvi.txt')
PAR_total_dict = datafile(outpath,lat,lon,'PAR_total.txt')
PAR_diffuse_dict = datafile(outpath,lat,lon,'PAR_diffuse.txt')
CF_dict =  datafile(outpath,lat,lon,'cloudfrac.txt')
toa_allsky_dict =  datafile(outpath,lat,lon,'toa_allsky.txt')
toa_par_tot_allsky_dict = datafile(outpath,lat,lon,'toa_par_tot_allsky.txt')

#extracting the monthly data and corresponding times 
AOD, AOD_time = AOD_dict['AOD Measurements'], AOD_dict['Times']
PRECIP, PRECIP_time = PRECIP_dict['Precipitation'], PRECIP_dict['Times']
NDVI, NDVI_time = NDVI_dict['NDVIs'], NDVI_dict['Times']
PAR_total, PAR_total_time = PAR_total_dict['PARs'],PAR_total_dict['Times']
PAR_diffuse, PAR_diffuse_time = PAR_diffuse_dict['PARs'], PAR_diffuse_dict['Times']
CF, CF_time = CF_dict['cloud_fracs'] ,CF_dict['Times']
toa_allsky, toa_allsky_time = toa_allsky_dict['PARs'] ,toa_allsky_dict['Times']
daily_solar_toa = toa_allsky_dict['Daily Solar Radiation']
toa_par_tot_allsky, toa_par_tot_allsky_time = toa_par_tot_allsky_dict['PARs'], toa_par_tot_allsky_dict['Times']

#standardised PAR anomalies
PAR_total_standardised_anoms = PAR_total_dict['Standardised Anomalies']
PAR_diffuse_standardised_anoms = PAR_diffuse_dict['Standardised Anomalies']

PAR_direct = []
Direct_std_anomalies = []

#calculating PAR direct from PAR total - PAR diffuse 

for i in range(len(PAR_diffuse)):
    PAR_direct.append(PAR_total[i] - PAR_diffuse[i])
    Direct_std_anomalies.append(PAR_total_standardised_anoms[i] - PAR_diffuse_standardised_anoms[i])

PAR_direct = np.asarray(PAR_direct)
Direct_std_anomalies = np.asarray(Direct_std_anomalies)

gamma = (daily_solar_toa/toa_allsky)

PAR_diff_corr = PAR_diffuse*gamma
PAR_tot_corr = PAR_total*gamma
PAR_direct_corr = PAR_direct*gamma

#using the scatterplot function in pltplots.py script 
# scatterplot(x_data, y_data, xlabel, ylabel,title)
scatterplot(toa_allsky_time,daily_solar_toa,'Time (month)','Daily Solar toa (W/$m^2$', 'Daily Solar toa 10S 64W')

scatterplot(PAR_total_time,PAR_tot_corr, 'Time (month)','PAR total corrected', 'Total PAR - Monthly average 10S 64W')

scatterplot(toa_allsky_time,toa_allsky,'Time (month)','toa allsky', ' TOA_Allsky Monthly 10S 64W' )

scatterplot(PAR_total_time,toa_par_tot_allsky,'Time (month)','toa PAR total allsky', 'PAR toa total allsky 10S 64W') 

#using pltplot.py lineofbestfit function to produce a scatter plot with a line of best fit and saving image to a .pdf file
#lineofbestfit(xdata,ydata, xlabel, ylabel, savefigure_path, figure suffix)
lineofbestfit(CF,PAR_diffuse, 'Cloud Fraction', 'PAR Diffuse (W/m^2)', figpath,lat+lon+ 'par_diff_CF.pdf')
lineofbestfit(CF,PAR_total, 'Cloud Fraction', 'PAR Total (W/m^2)', figpath, lat+lon+'par_diff_CF.pdf')
lineofbestfit(CF,PAR_direct, 'Cloud Fraction', 'PAR Direct (W/m^2)', figpath, lat+lon+'par_diff_CF.pdf')

# one plot of CF against the three datasets: PAR_diffuse, PAR_total and PAR_direct
triplotsharex(CF,CF,CF,PAR_diffuse,PAR_total,PAR_direct,'Cloud Fraction','PAR_diff', 'PAR_tot', 'PAR_direct', figpath,lat+lon+"CF_par_corr.pdf",'PAR vs Cloud Fraction ' +lat+ "N " +lon + "E")



#There are fewer AOD data points than for PAR over the same time period
#xtime_equals_ytime is a function that filters through the x and y time lists, locates when the times are equal
#and creates a numpy array of index values when xtime == ytime
#index array can be used to extract the data points in the longer list for the times when both datasets have readings

PAR_diff_at_AOD =  PAR_diffuse[xtime_equals_ytime(AOD_time,PAR_diffuse_time)]
PAR_total_at_AOD = PAR_total[xtime_equals_ytime(AOD_time,PAR_total_time)]
PAR_direct_at_AOD = PAR_direct[xtime_equals_ytime(AOD_time,PAR_total_time)]


# three plot of AOD against the datasets: PAR_diffuse, PAR_total and PAR_direct - NOT CORRECTED
#x axis is shared
#triploshare(x1,x2,x3,y1,y2,y3,xlabel, y1label, y2label, y3label, figurepath, figure suffix, figure title )
triplotsharex(AOD,AOD,AOD,PAR_diff_at_AOD,PAR_total_at_AOD,PAR_direct_at_AOD, 'AOD', 'PAR_diff / W/$m^2$', 'PAR_tot / W/$m^2$', 'PAR_direct / W/$m^2$',figpath,lat+lon+ "AOD_par.pdf", "AOD vs PAR " +lat+ "N " +lon + "E") 


# three plots of CF against datasets: PAR_diffuse, PAR_total and PAR_direct - CORRECTED
triplotsharex(CF,CF,CF,PAR_diff_corr,PAR_tot_corr,PAR_direct_corr, 'Cloud Fraction', 'PAR_diff_corr / W/$m^2$','PAR_tot_corr / W/$m^2$','PAR_direct_corr / W/$m^2$',figpath, lat+lon+ "CF_par_corr.pdf", "Cloud Fraction vs PAR_Corrected " +lat+ "N " +lon + "E" )


PAR_diff_at_AOD_corr = PAR_diff_corr[xtime_equals_ytime(AOD_time,PAR_diffuse_time)]
PAR_tot_at_AOD_corr = PAR_tot_corr[xtime_equals_ytime(AOD_time,PAR_total_time)]
PAR_direct_at_AOD_corr =  PAR_direct_corr[xtime_equals_ytime(AOD_time,PAR_total_time)]

# Three plos AOD against the three datasets: PAR_diffuse, PAR_total and PAR_direct - CORRECTED
triplotsharex(AOD,AOD,AOD,PAR_diff_at_AOD_corr ,PAR_tot_at_AOD_corr,PAR_direct_at_AOD_corr, 'AOD', 'PAR_diff / W/$m^2$', 'PAR_tot / W/$m^2$', 'PAR_direct / W/$m^2$',figpath,lat+lon+"AOD_par_corected.pdf", "AOD vs PAR Corrected "+lat+ "N " +lon + "E")


#3x2 plot of the monthly data readings against time
#line of best fit 
#yearly average also plotted
gridplot(AOD_time, CF_time, NDVI_time, PAR_diffuse_time, PAR_total_time, PRECIP_time,AOD, CF, NDVI, PAR_diffuse, PAR_total, PRECIP,figpath,lat+lon+"monthly_data_v_time.pdf", "Monthly Data (2002 - 2012) "+lat+ "N " + lon + "E") 


#alternative way of producing one plot with a shared x axis 

fig4 = plt.figure(figsize=(4,4))
    
plt.plot(CF,PAR_diffuse,'o',label = 'PAR diffuse')
plt.xlabel('Cloud Fraction')
plt.ylabel('PAR (W/$m^2$)')
bestfitline(CF,PAR_diffuse)
plt.plot(CF,PAR_total,'o', label = 'PAR total')
bestfitline(CF,PAR_total)
plt.plot(CF,PAR_direct,'o', label = 'PAR direct')
bestfitline(CF,PAR_direct)
plt.show()
fig4.savefig(figpath + lat+lon+"CF_par_controlregion.pdf" , bbox_inches='tight')


#top of atmosphere allsky plot against cloud fraction for january
toa_at_AOD = toarray(toa_allsky)[xtime_equals_ytime(toa_allsky_time,AOD_time)]


#extracting only the January points from PAR using the decimal place 
#e.g. TIME = YEAR + (MONTH-0.5/12) so 
#therefore JAN = 1-0.5/12 = 0.042
#%1.0 extracts the decimal from a list of floats
  
Jan_PAR_ID = np.where((toa_allsky_time %1.0 > 0) & (toa_allsky_time %1.0 < 0.1))
Jan_CF_ID = np.where((CF_time %1.0 > 0) & (CF_time %1.0 < 0.1))


plt.plot(toa_allsky_time[Jan_PAR_ID],toa_allsky[Jan_PAR_ID], 'o')
plt.xlabel('Year of Jan Retrival')
plt.ylabel('PAR toa all sky (W/m^2)')
plt.show()

plt.plot(CF[Jan_CF_ID],toa_allsky[Jan_PAR_ID], 'o')
plt.xlabel('Cloud Fraction Jan')
plt.ylabel('PAR toa allsky Jan')
plt.show()

plt.plot( CF, toa_allsky,'o', label = 'cloud fraction')
plt.xlabel('Cloud Fraction')
plt.ylabel('toa_allsky')
plt.legend()
plt.show()


April_IDaod = np.where((AOD_time %1.0 > 0.21) & (AOD_time %1.0 <0.3))
AOD_April_times = np.array(AOD_time[April_IDaod])
AOD_April = np.array(AOD[April_IDaod])


April_IDpar_times = np.where((PAR_diffuse_time%1.0 > 0.21) & (PAR_diffuse_time %1.0 <0.3))

PAR_April_times = np.array(PAR_diffuse_time[April_IDpar_times]) 
PAR_April_diff = np.array(PAR_diffuse[April_IDpar_times])
PAR_April_tot = np.array(PAR_total[April_IDpar_times])
PAR_April_direct = np.array(PAR_direct[April_IDpar_times])

diffuse = PAR_April_diff[xtime_equals_ytime(PAR_April_times,AOD_April_times)]
total = PAR_April_tot[xtime_equals_ytime(PAR_April_times,AOD_April_times)]
direct = PAR_April_direct[xtime_equals_ytime(PAR_April_times,AOD_April_times)]


lineofbestfit(AOD_April,direct,'AOD', 'PAR direct April (W/m^2)',figpath, lat+lon+'aod_par_direct_april.pdf')
lineofbestfit(AOD_April, total, 'AOD','PAR total April (W/m^2)', figpath, lat+lon+ 'aod_par_total_april.pdf')
lineofbestfit(AOD_April,diffuse,'AOD','PAR diffuse April (W/m^2)',figpath, lat+ lon+ 'aod_par_diff_april.pdf')




#--------------------------------
#separating data into rainy season and dry season before plotting
#rainy season in the amazon occurs between December and May
#dry season is from June to November

year_list_AOD = linspace(2002,2012,11)
months = linspace(1,12,12)
year_list_PRECIP = linspace(2000,2018,19)

aoddry = toarray(AOD_dict['dry anomalies'])
aodwet = toarray(AOD_dict['rainy anomalies'])
aodwet_t = toarray(AOD_dict['rainy times'])
aoddry_t = toarray(AOD_dict['dry times'])

cloud_fracdry = toarray(CF_dict['dry anomalies'])
cloud_fracwet = toarray(CF_dict['rainy anomalies'])
cloud_fracwet_t = toarray(CF_dict['rainy times'])
cloud_fracdry_t = toarray(CF_dict['dry times'])


NDVIdry = toarray(NDVI_dict['dry anomalies'])
NDVIwet = toarray(NDVI_dict['rainy anomalies'])
NDVIwet_t = toarray(NDVI_dict['rainy times'])
NDVIdry_t = toarray(NDVI_dict['dry times'])

NDVI_yearly_anom_dry = toarray(NDVI_dict['nine year anoms wet'])    
NDVI_yearly_anom_wet = toarray(NDVI_dict['nine year anoms dry']) 

cloud_frac_yearly_anom_dry = toarray(CF_dict['nine year anoms wet'])    
cloud_frac_yearly_anom_wet = toarray(CF_dict['nine year anoms dry']) 

aod_yearly_anom_dry = toarray(AOD_dict['nine year anoms wet'])    
aod_yearly_anom_wet = toarray(AOD_dict['nine year anoms dry']) 

lineofbestfit(aod_yearly_anom_dry,cloud_frac_yearly_anom_dry, 'aod dry season yearly anomaly', 'Cloud fraction dry season anomalies',figpath, lat+lon+'aod_cloud_frac_rs.pdf')
lineofbestfit(aod_yearly_anom_wet,cloud_frac_yearly_anom_wet, 'aod rainy season yearly anomaly', 'Cloud fraction rainy season anomalies',figpath,  lat+lon+'aod_cloud_frac_rs.pdf')
lineofbestfit(NDVI_yearly_anom_wet,cloud_frac_yearly_anom_wet, 'NDVI rainy season yearly anomaly', 'Cloud fraction rainy season anomalies',figpath,  lat+lon+'ndvi_cloud_frac_rs.pdf')
lineofbestfit(NDVI_yearly_anom_dry,cloud_frac_yearly_anom_dry, 'NDVI dry season yearly anomaly', 'Cloud fraction dry season anomalies',figpath,  lat+lon+'ndvi_cloud_frac_ds.pdf')
lineofbestfit(NDVIwet,cloud_fracwet, 'NDVI rainy season anomalies', 'Cloud fraction rainy season anomalies',figpath, lat+lon+'ndvi_cloud_frac_rs.pdf')


cloud_frac_at_aod_dry = cloud_fracdry[xtime_equals_ytime(aoddry_t,cloud_fracdry_t)]
cloud_frac_at_aod_wet = cloud_fracwet[xtime_equals_ytime(aodwet_t,cloud_fracwet_t)]
lineofbestfit(aodwet,cloud_frac_at_aod_wet, 'aod rainy season anomalies', 'Cloud fraction rainy season anomalies',figpath, lat+lon+ 'ndvi_cloud_frac_rs.pdf')


#simple scatter plots of two data sets of equal length

plt.plot(cloud_frac_at_aod_dry, aoddry, 'o')
plt.xlabel('CF dry montly data')
plt.ylabel('AOD wet monthly data')
plt.show()
plt.plot(cloud_frac_at_aod_wet, aodwet, 'o')
plt.xlabel('CF wet monthly data')
plt.ylabel('AOD wet monthly data')

plt.show()


#function to convert a list of floats into their log equivalent
def tolog(mylist):
    logs = []
    for i in range(len(mylist)):
        log = math.log10((mylist)[i])
        logs.append(log)
    return logs




x1 = toarray(AOD_dict['Times'])
x2 = toarray(NDVI_dict['Times'])
x3 = toarray(PRECIP_dict['Times'])
x4a = toarray(PAR_total_dict['Times'])
x4b = toarray(PAR_diffuse_dict['Times'])
x5 = toarray(CF_dict['Times'])

y1 = toarray(AOD_dict['Anomalies'])
y2 = toarray(NDVI_dict['Anomalies'])
y3 = toarray(PRECIP_dict['Anomalies'])
y4a = toarray(PAR_total_dict['Anomalies'])
y4b = toarray(PAR_diffuse_dict['Anomalies'])
y5 = toarray(CF_dict['Anomalies'])


p = gridplot(x1,x2,x3,x4a,x4b,x5,y1,y2,y3,y4a,y4b,y5,figpath,lat+lon+"monthlyanomalies.pdf","Time series plots of anomalies for " + lat + "N " + lon + "E")

#--------------------------------------------------

#precipitation vs cloud fraction
#precipitation data is larger than cf data

precip_at_cf = PRECIP[xtime_equals_ytime(x3,x5)]
lineofbestfit(precip_at_cf, CF, 'recipitation / mm/hr ', 'Cloud Fraction',figpath, lat+lon +'cfvprecip.pdf')

#NDVI vs PAR
ndvi_at_PAR = NDVI[xtime_equals_ytime(x2,x4b)]
lineofbestfit(ndvi_at_PAR,PAR_diffuse, 'NDVI','PAR_diffuse / W/$m^2$',figpath, lat+lon+'ndvi_PAR_d.pdf')

#NDVI vs CF
ndvi_at_cf = NDVI[xtime_equals_ytime(x2,x5)]
lineofbestfit(ndvi_at_cf,CF,'NDVI','Cloud Fraction', figpath, lat+lon+'ndvi_cd.pdf')

#AOD vs Precipitation
precip_at_AOD = PRECIP[xtime_equals_ytime(x1,x3)]
lineofbestfit(precip_at_AOD,AOD, 'Precipitation / mm/hr ', 'AOD', figpath, lat+lon+'precip_aod.pdf')



