#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
NAME;
Script to process AOD, Precipitation and NDVI data for a coodinate grid using sub processing scripts 
and use matplotlib to create time series plots

PURPOSE;
This program inputs the lat lon cooridnate range under consideration, and outputs four plots:
1. Variable(AOD,NDVI Precipiation or cloud_frac) vs month
2. Number of variable retrievals over the 17 year period vs month
3. Average monthly variable reading vs month
4. Variable anomalies vs year 

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

lon = -63.5 
lat = -10.5 


#lat_bnds, lon_bnds = [-12, -11], [-64, -63]

lat_bnds, lon_bnds = [-11,-9], [-65,-63]

#lat and lon bounds for region which hasn't undergone deforestation
lat_bnds_controlregion, lon_bnds_controlregion = [-9,-7], [-67,-65]


#processing data for the region which has not undergone deforestation in the Amazon, North of Rondonia and Mato Grosso
Aerosol_data_controlregion = process_aerosol(lon_bnds_controlregion,lat_bnds_controlregion,outpath,aerosol_path,'aodfiles_control.json', 'aodlonlat_control.csv', 'AOD_no_nan_control.json')

Precipitation_data_controlregion = process_precip(lon_bnds_controlregion,lat_bnds_controlregion,outpath,precip_path,"precip_data_control.txt", "mean_precip_data_control.csv")

cloud_frac_data_controlregion = process_cloud_frac(lat_bnds_controlregion,lon_bnds_controlregion,outpath, cloud_frac_path, "v2.0.nc", "cf_file_info_controlregion.csv", "cf_data_controlregion.json")

NDVI_data_controlregion = process_ndvi(lon_bnds_controlregion,lat_bnds_controlregion,outpath,ndvi_path,".hdf","NDVI_files_controlregion.json", "mean_NDVI_data_controlregion.csv")

PAR_Total_data_controlregion = process_PAR(outpath,PAR_path,lat_bnds_controlregion,lon_bnds_controlregion,"boa_par_tot_allsky","fv2.1.nc","PAR_tot_file_info_controlregion.csv", "PAR_tot_data_controlregion.json") #base of atmosphere total allsky composite w/m^2
PAR_Diffuse_data_controlregion = process_PAR(outpath,PAR_path,lat_bnds_controlregion,lon_bnds_controlregion,"boa_par_dif_allsky", "fv2.1.nc","PAR_diff_file_info_controlregion.csv", "PAR_diff_data_controlregion.json" )#base of atmosphere diffuse allsky composite w/m^2

toa_allsky_data_controlregion = process_PAR(outpath,PAR_path,lat_bnds_controlregion,lon_bnds_controlregion, "toa_swdn_allsky","fv2.1.nc", "toa_allsky_vals_controlregion.csv", "toa_attributes.json")
toa_par_tot_allsky_data_controlregion = process_PAR(outpath,PAR_path,lat_bnds_controlregion, lon_bnds_controlregion, "toa_par_tot_allsky","fv2.1.nc", "toa_par_tot_allsky_vals_controlregion.csv", "toa_par_tot_attributes_controlregion")



#processing data for Rondonia where rapid deforestation has occured since ~ 1990 
Aerosol_data = process_aerosol(lon_bnds,lat_bnds,outpath,aerosol_path,'aodfiles.json','aodlonlat.csv','AOD_no_nan.json')

Precipitation_data = process_precip(lon_bnds,lat_bnds,outpath,precip_path,'precipitation.txt', 'precip_means.csv') 

NDVI_data = process_ndvi(lon_bnds,lat_bnds,outpath,ndvi_path,".hdf","NDVI_files.json", "mean_ndvi_out.csv")

PAR_Total_data = process_PAR(outpath,PAR_path,lat_bnds,lon_bnds,"boa_par_tot_allsky","fv2.1.nc","PAR_total_vals.csv","par_total_attributes.json") #base of atmosphere total allsky composite w/m^2
PAR_Diffuse_data = process_PAR(outpath,PAR_path,lat_bnds,lon_bnds,"boa_par_dif_allsky", "fv2.1.nc", "PAR_diffuse_vals.csv","par_diffuse_attributes.json" )#base of atmosphere diffuse allsky composite w/m^2

cloud_frac_data = process_cloud_frac(lat_bnds,lon_bnds,outpath, cloud_frac_path, "v2.0.nc","cloud_frac_vals.csv", "cloud_frac_attributes.json")

toa_allsky_data = process_PAR(outpath,PAR_path,lat_bnds, lon_bnds, "toa_swdn_allsky","fv2.1.nc", "toa_allsky_vals.csv", "toa_attributes.json")
toa_par_tot_allsky_data = process_PAR(outpath,PAR_path,lat_bnds, lon_bnds, "toa_par_tot_allsky","fv2.1.nc", "toa_par_tot_allsky_vals.csv", "toa_par_tot_allsky_attributes.json" )


#monthly data readings and time of reading for PAR, CF,AOD, toa_allsky and toa_pa_tot_allsky for Rondonia
PAR_diff, PAR_diff_times = PAR_Diffuse_data['PARs'],PAR_Diffuse_data['Times']
PAR_tot,PAR_tot_times = PAR_Total_data['PARs'],PAR_Diffuse_data['Times']

CF,CF_times = cloud_frac_data['cloud_fracs'],cloud_frac_data['Times']

AOD, Aerosol_times = Aerosol_data['AOD Measurements'], Aerosol_data['Times']

toa_allsky,toa_allsky_times = toa_allsky_data['PARs'],toa_allsky_data['Times']

NDVI,NDVI_times = NDVI_data['NDVIs'], NDVI_data['Times']
#Theoretical solar flux
daily_solar_toa = toa_allsky_data['Daily Solar Radiation']
toa_par_tot_allsky = toa_par_tot_allsky_data['PARs']




#monthly data readings and time of reading for the control region north of Rondonia

PAR_diff_CONTROL,PAR_diff_CONTROL_time = PAR_Diffuse_data_controlregion['PARs'],PAR_Diffuse_data_controlregion['Times']
PAR_tot_CONTROL,PAR_tot_CONTROL_time = PAR_Total_data_controlregion['PARs'],PAR_Total_data_controlregion['Times']

NDVI_CONTROL,NDVI_CONTROL_time = NDVI_data_controlregion['NDVIs'],NDVI_data_controlregion['Times']
AOD_CONTROL, AOD_CONTROL_time = Aerosol_data_controlregion['AOD Measurements'],Aerosol_data_controlregion['Times']

Precip_CONTROL,Precip_CONTROL_time = Precipitation_data_controlregion['Precipitation'], Precipitation_data_controlregion['Times']

CF_CONTROL,CF_CONTROL_time = cloud_frac_data_controlregion['cloud_fracs'], cloud_frac_data_controlregion['Times']

toa_allsky_CONTROL,toa_allsky_CONTROL_time = toa_allsky_data_controlregion['PARs'],toa_allsky_data_controlregion['Times']

toa_par_toa_allsky_CONTROL,toa_par_tot_allsky_CONTROL_time = toa_par_tot_allsky_data_controlregion['PARs'],toa_par_tot_allsky_data_controlregion['Times']

daily_solar_toa_CONTROL = toa_allsky_data_controlregion['Daily Solar Radiation']

#standardised PAR anomalies

PAR_tot_standardised_anoms = PAR_Total_data['Standardised Anomalies']
PAR_diffuse_standardised_anoms = PAR_Diffuse_data['Standardised Anomalies']

PAR_direct = []
Direct_std_anomalies = []

#calculating PAR direct from PAR total - PAR diffuse 
for i in range(len(PAR_diff)):
    PAR_direct.append(PAR_tot[i] - PAR_diff[i])
    Direct_std_anomalies.append(PAR_tot_standardised_anoms[i] - PAR_diffuse_standardised_anoms[i])

PAR_direct = np.asarray(PAR_direct)
Direct_std_anomalies = np.asarray(Direct_std_anomalies)

PAR_direct_CONTROL = []
for i in range(len(PAR_diff_CONTROL)):
    PAR_direct_CONTROL.append(PAR_tot_CONTROL[i] - PAR_diff_CONTROL[i])

PAR_direct_CONTROL = np.asarray(PAR_direct_CONTROL)

#applying the correction which arises from daily solar fluctuations to PAR measurements

gamma = (daily_solar_toa/toa_allsky)

PAR_diff_corr = PAR_diff*gamma
PAR_tot_corr = PAR_tot*gamma
PAR_direct_corr = PAR_direct*gamma

#applying corection to control region location

gamma2 = (daily_solar_toa_CONTROL/toa_allsky_CONTROL)
PAR_diff_corr_CONTROL = PAR_diff_CONTROL*gamma2
PAR_tot_corr_CONTROL = PAR_tot_CONTROL*gamma2
PAR_direct_corr_CONTROL = PAR_direct_CONTROL*gamma2



scatterplot(toa_allsky_times,daily_solar_toa,'Time (month)','Daily Solar toa (W/$m^2$', 'Daily Solar toa 10S 64W')

scatterplot(PAR_tot_times,PAR_tot_corr, 'Time (month)','PAR total corrected', 'Total PAR - Monthly average 10S 64W')

scatterplot(toa_allsky_times,toa_allsky,'Time (month)','toa allsky', ' TOA_Allsky Monthly 10S 64W' )

scatterplot(PAR_tot_times,toa_par_tot_allsky,'Time (month)','toa PAR total allsky', 'PAR toa total allsky 10S 64W') 

scatterplot(NDVI_CONTROL_time,NDVI_CONTROL,'Time (month)', 'NDVI', 'NDVI Monthly 8S 66W')

scatterplot(NDVI_times,NDVI,'Time(month)', 'NDVI', 'NDVI Monthly 10S 64W')


gridplot(AOD_CONTROL_time, CF_CONTROL_time, NDVI_CONTROL_time, PAR_diff_CONTROL_time, PAR_tot_CONTROL_time, Precip_CONTROL_time,AOD_CONTROL, CF_CONTROL, NDVI_CONTROL, PAR_diff_CONTROL, PAR_tot_CONTROL, Precip_CONTROL,figpath,"monthly_data_v_time_controlregion.pdf", "Monthly Data (2002 - 2012) 8S, 66W") 



#calculating standardised anomalies for PAR (diffuse,total and direct) and cloud fraction
#plotting PAR (corrected and non-corrected for daily variations in solar flux at toa) against cloud fraction for Rondonia

lineofbestfit(CF,PAR_diff, 'Cloud Fraction', 'PAR Diffuse (W/m^2)', figpath, 'par_diff_CF_r.pdf')
lineofbestfit(CF,PAR_tot, 'Cloud Fraction', 'PAR Total (W/m^2)', figpath, 'par_diff_CF_r.pdf')
lineofbestfit(CF,PAR_direct, 'Cloud Fraction', 'PAR Direct (W/m^2)', figpath, 'par_diff_CF_r.pdf')


# one plot of CF against the three datasets: PAR_diffuse, PAR_total and PAR_direct - CORRECTED
triplotsharex(CF,CF,CF,PAR_diff,PAR_tot,PAR_direct,'Cloud Fraction','PAR_diff', 'PAR_tot', 'PAR_direct', figpath,"CF_par_corr_r.pdf",'PAR vs Cloud Fraction 10S 64W')

# one plot of CF against the three datasets: PAR_diffuse, PAR_total and PAR_direct - CORRECTED

triplotsharex(CF,CF,CF,PAR_diff_corr,PAR_tot_corr,PAR_direct_corr, 'Cloud Fraction', 'PAR_diff_corr / W/$m^2$','PAR_tot_corr / W/$m^2$','PAR_direct_corr / W/$m^2$',figpath, "CF_par_corr_r.pdf", "Cloud Fraction vs PAR_Corrected 10S 64W")

# one plot of AOD against the three datasets: PAR_diffuse, PAR_total and PAR_direct - NOT CORRECTED 

triplotsharex(AOD,AOD,AOD,PAR_diff[xtime_equals_ytime(Aerosol_times,PAR_diff_times)],PAR_tot[xtime_equals_ytime(Aerosol_times,PAR_tot_times)],PAR_direct[xtime_equals_ytime(Aerosol_times,PAR_tot_times)], 'AOD', 'PAR_diff / W/$m^2$', 'PAR_tot / W/$m^2$', 'PAR_direct / W/$m^2$',figpath,"AOD_par_r.pdf", "AOD vs PAR") 


# one plot of AOD against the three datasets: PAR_diffuse, PAR_total and PAR_direct - CORRECTED
triplotsharex(AOD,AOD,AOD,PAR_diff_corr[xtime_equals_ytime(Aerosol_times,PAR_diff_times)],PAR_tot_corr[xtime_equals_ytime(Aerosol_times,PAR_tot_times)],PAR_direct_corr[xtime_equals_ytime(Aerosol_times,PAR_tot_times)], 'AOD', 'PAR_diff / W/$m^2$', 'PAR_tot / W/$m^2$', 'PAR_direct / W/$m^2$',figpath,"AOD_par_corected_r.pdf", "AOD vs PAR Corrected")


#plotting PAR (corrected and non-corrected for daily variations in solar flux at toa) against cloud fraction
#for a region directly North of Rondonia which has not undergone deforesation


fig4 = plt.figure(figsize=(4,4))
    
plt.plot(CF,PAR_diff,'o',label = 'PAR diffuse')
plt.xlabel('Cloud Fraction')
plt.ylabel('PAR (W/$m^2$)')
bestfitline(CF_CONTROL,PAR_diff_CONTROL)
plt.plot(CF_CONTROL,PAR_tot_CONTROL,'o', label = 'PAR total')
bestfitline(CF_CONTROL,PAR_tot_CONTROL)
plt.plot(CF_CONTROL,PAR_direct_CONTROL,'o', label = 'PAR direct')
bestfitline(CF_CONTROL,PAR_direct_CONTROL)
plt.show()
fig4.savefig(figpath + "CF_par_controlregion.pdf", bbox_inches='tight')

fig5 = plt.figure(figsize=(4,4))
plt.plot(CF_CONTROL,PAR_diff_corr_CONTROL,'o',label = 'PAR diffuse')
plt.xlabel('Cloud Fraction')
plt.ylabel('PAR_Corrected (W/$m^2$)')
bestfitline(CF_CONTROL,PAR_diff_corr_CONTROL)
plt.plot(CF_CONTROL,PAR_tot_corr_CONTROL,'o', label = 'PAR total')
bestfitline(CF_CONTROL,PAR_tot_corr_CONTROL)
plt.plot(CF_CONTROL,PAR_direct_corr_CONTROL,'o', label = 'PAR direct')
bestfitline(CF_CONTROL,PAR_direct_corr_CONTROL)
plt.show()
fig5.savefig(figpath + "CF_par_corr_controlregion.pdf", bbox_inches='tight')


fig6 = plt.figure(figsize=(4,4))
plt.plot(AOD_CONTROL,PAR_diff_CONTROL[xtime_equals_ytime(AOD_CONTROL_time,PAR_diff_CONTROL_time)],'o',label = 'PAR diffuse')
plt.xlabel('AOD')
plt.ylabel('PAR(W/$m^2$)')
bestfitline(AOD_CONTROL,PAR_diff_CONTROL[xtime_equals_ytime(AOD_CONTROL_time,PAR_diff_CONTROL_time)])
plt.plot(AOD_CONTROL,PAR_tot_CONTROL[xtime_equals_ytime(AOD_CONTROL_time,PAR_tot_CONTROL_time)],'o', label = 'PAR total')
bestfitline(AOD_CONTROL,PAR_tot_CONTROL[xtime_equals_ytime(AOD_CONTROL_time,PAR_tot_CONTROL_time)])
plt.plot(AOD_CONTROL,PAR_direct_CONTROL[xtime_equals_ytime(AOD_CONTROL_time,PAR_tot_CONTROL_time)],'o', label = 'PAR direct')
bestfitline(AOD_CONTROL,PAR_direct_CONTROL[xtime_equals_ytime(AOD_CONTROL_time,PAR_tot_CONTROL_time)])
plt.show()
fig6.savefig(figpath + "AOD_par_r.pdf", bbox_inches='tight')

fig7 = plt.figure(figsize=(4,4))
plt.plot(AOD_CONTROL,PAR_diff_corr_CONTROL[xtime_equals_ytime(AOD_CONTROL_time,PAR_diff_CONTROL_time)],'o',label = 'PAR diffuse')
plt.xlabel('AOD')
plt.ylabel('PAR_Corrected (W/$m^2$)')
bestfitline(AOD_CONTROL,PAR_diff_corr_CONTROL[xtime_equals_ytime(AOD_CONTROL_time,PAR_diff_CONTROL_time)])
plt.plot(AOD_CONTROL,PAR_tot_corr_CONTROL[xtime_equals_ytime(AOD_CONTROL_time,PAR_tot_CONTROL_time)],'o', label = 'PAR total')
bestfitline(AOD_CONTROL,PAR_tot_corr_CONTROL[xtime_equals_ytime(AOD_CONTROL_time,PAR_tot_CONTROL_time)])
plt.plot(AOD_CONTROL,PAR_direct_corr_CONTROL[xtime_equals_ytime(AOD_CONTROL_time,PAR_tot_CONTROL_time)],'o', label = 'PAR direct')
bestfitline(AOD_CONTROL,PAR_direct_corr_CONTROL[xtime_equals_ytime(AOD_CONTROL_time,PAR_tot_CONTROL_time)])
plt.show()
fig7.savefig(figpath + "AOD_par_corr_r.pdf", bbox_inches='tight')

sys.exit()
#top of atmosphere allsky plot against cloud fraction for january
toa_at_AOD = toarray(toa_allsky)[xtime_equals_ytime(toa_allsky_times,Aerosol_times)]

Jan_PAR_ID = np.where((toa_allsky_times %1.0 > 0) & (toa_allsky_times %1.0 < 0.1))
Jan_CF_ID = np.where((CF_times %1.0 > 0) & (CF_times %1.0 < 0.1))

plt.plot(toa_allsky_times[Jan_PAR_ID],toa_allsky[Jan_PAR_ID], 'o')
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


#plotting PAR against AOD

lineofbestfit(Direct_std_anomalies[xtime_equals_ytime(PAR_diff_times,Aerosol_times)],Aerosol_data['Standardised Anomalies'],'AOD', 'PAR direct (W/m^2)',figpath, 'aod_par_direct_std_r.pdf')
lineofbestfit(np.asarray(PAR_tot_standardised_anoms)[xtime_equals_ytime(PAR_diff_times,Aerosol_times)],Aerosol_data['Standardised Anomalies'],'AOD', 'PAR total (W/m^2)',figpath, 'aod_par_tot_std_r.pdf')
lineofbestfit(np.asarray(PAR_diffuse_standardised_anoms)[xtime_equals_ytime(PAR_diff_times,Aerosol_times)],Aerosol_data['Standardised Anomalies'],'AOD', 'PAR diffuse (W/m^2)',figpath, 'aod_par_diffuse_std_r.pdf')

April_IDaod = np.where((Aerosol_times %1.0 > 0.21) & (Aerosol_times %1.0 <0.3))


AOD_April_times = np.array(Aerosol_times[April_IDaod])
AOD_April = np.array(AOD[April_IDaod])


April_IDpar_times = np.where((PAR_diff_times%1.0 > 0.21) & (PAR_diff_times %1.0 <0.3))
PAR_April_times = np.array(PAR_diff_times[April_IDpar_times]) 
PAR_April_diff = np.array(PAR_diff[April_IDpar_times])
PAR_April_tot = np.array(PAR_tot[April_IDpar_times])
PAR_April_direct = np.array(PAR_direct[April_IDpar_times])

diffuse = PAR_April_diff[xtime_equals_ytime(PAR_April_times,AOD_April_times)]
total = PAR_April_tot[xtime_equals_ytime(PAR_April_times,AOD_April_times)]
direct = PAR_April_direct[xtime_equals_ytime(PAR_April_times,AOD_April_times)]


lineofbestfit(AOD_April,direct,'AOD', 'PAR direct April (W/m^2)',figpath, 'aod_par_direct_april_r.pdf')
lineofbestfit(AOD_April, total, 'AOD','PAR total April (W/m^2)', figpath, 'aod_par_total_april_r.pdf')
lineofbestfit(AOD_April,diffuse,'AOD','PAR diffuse April (W/m^2)',figpath, 'aod_par_diff_april_r.pdf')


#lineofbestfit(AOD,PARs_direct,  'AOD', 'PAR direct (W/m^2)',figpath, 'aod_par_direct_monthly_r.pdf')
#lineofbestfit(AOD, PARs_total, 'AOD','PAR total (W/m^2)', figpath, 'aod_par_total_monthly_r.pdf')
#lineofbestfit(AOD,PARs_diffuse,   'AOD','PAR diffuse(W/m^2)',figpath, 'aod_par_direct_monthly_r.pdf')


#lineofbestfit(Aerosol_data,PARs_direct,  'AOD', 'PAR direct (W/m^2)',figpath, 'aod_par_direct_monthly_r.pdf')

        

#--------------------------------
#plotting
#1) separating data into rainy season and dry season before plotting

year_list_AOD = linspace(2002,2012,11)
months = linspace(1,12,12)
year_list_PRECIP = linspace(2000,2018,19)

aoddry = toarray(Aerosol_data['dry anomalies'])
aodwet = toarray(Aerosol_data['rainy anomalies'])
aodwet_t = toarray(Aerosol_data['rainy times'])
aoddry_t = toarray(Aerosol_data['dry times'])

cloud_fracdry = toarray(cloud_frac_data['dry anomalies'])
cloud_fracwet = toarray(cloud_frac_data['rainy anomalies'])
cloud_fracwet_t = toarray(cloud_frac_data['rainy times'])
cloud_fracdry_t = toarray(cloud_frac_data['dry times'])


NDVIdry = toarray(NDVI_data['dry anomalies'])
NDVIwet = toarray(NDVI_data['rainy anomalies'])
NDVIwet_t = toarray(NDVI_data['rainy times'])
NDVIdry_t = toarray(NDVI_data['dry times'])

NDVI_yearly_anom_dry = toarray(NDVI_data['nine year anoms wet'])    
NDVI_yearly_anom_wet = toarray(NDVI_data['nine year anoms dry']) 

cloud_frac_yearly_anom_dry = toarray(cloud_frac_data['nine year anoms wet'])    
cloud_frac_yearly_anom_wet = toarray(cloud_frac_data['nine year anoms dry']) 

aod_yearly_anom_dry = toarray(Aerosol_data['nine year anoms wet'])    
aod_yearly_anom_wet = toarray(Aerosol_data['nine year anoms dry']) 

lineofbestfit(aod_yearly_anom_dry,cloud_frac_yearly_anom_dry, 'aod dry season yearly anomaly', 'Cloud fraction dry season anomalies',figpath, 'aod_cloud_frac_rs.pdf')
lineofbestfit(aod_yearly_anom_wet,cloud_frac_yearly_anom_wet, 'aod rainy season yearly anomaly', 'Cloud fraction rainy season anomalies',figpath, 'aod_cloud_frac_rs.pdf')
lineofbestfit(NDVI_yearly_anom_wet,cloud_frac_yearly_anom_wet, 'NDVI rainy season yearly anomaly', 'Cloud fraction rainy season anomalies',figpath, 'ndvi_cloud_frac_rs.pdf')
lineofbestfit(NDVI_yearly_anom_dry,cloud_frac_yearly_anom_dry, 'NDVI dry season yearly anomaly', 'Cloud fraction dry season anomalies',figpath, 'ndvi_cloud_frac_ds.pdf')

plt.plot(NDVIdry,cloud_fracdry, 'o')
plt.xlabel('NDVI Dry')
plt.ylabel('CF Dry')
plt.show()

plt.plot(NDVIwet,cloud_fracwet, 'o')
plt.xlabel('NDVI wet')
plt.ylabel('CF Wet')
plt.show()

lineofbestfit(NDVIwet,cloud_fracwet, 'NDVI rainy season anomalies', 'Cloud fraction rainy season anomalies',figpath, 'ndvi_cloud_frac_rs.pdf')

#print len(cloud_fracdry), len(cloud_fracwet), len(aoddry), len(aodwet)
cloud_frac_at_aod_dry = cloud_fracdry[xtime_equals_ytime(aoddry_t,cloud_fracdry_t)]
cloud_frac_at_aod_wet = cloud_fracwet[xtime_equals_ytime(aodwet_t,cloud_fracwet_t)]
lineofbestfit(aodwet,cloud_frac_at_aod_wet, 'aod rainy season anomalies', 'Cloud fraction rainy season anomalies',figpath, 'ndvi_cloud_frac_rs.pdf')



plt.plot(cloud_frac_at_aod_dry, aoddry, 'o')
plt.xlabel('CF dry montly data')
plt.ylabel('AOD wet monthly data')
plt.show()
plt.plot(cloud_frac_at_aod_wet, aodwet, 'o')
plt.xlabel('CF wet monthly data')
plt.ylabel('AOD wet monthly data')

plt.show()


sys.exit()


#x1log = math.log10(Aerosol_data['Times'])

def tolog(mylist):
    logs = []
    for i in range(len(mylist)):
        log = math.log10((mylist)[i])
        logs.append(log)
    return logs

#print x1log, x1

aods = toarray(Aerosol_data['AOD Measurements'])
ndvis = toarray(tolog(NDVI_data['NDVIs']))
precips = toarray(tolog(Precipitation_data['Precipitation']))
cfs = toarray(cloud_frac_data['cloud_fracs'])
PAR_d = toarray(tolog(PAR_Diffuse_data['PARs']))


x1 = toarray(Aerosol_data['Times'])
x2 = toarray(NDVI_data['Times'])
x3 = toarray(Precipitation_data['Times'])
x4a = toarray(PAR_Total_data['Times'])
x4b = toarray(PAR_Diffuse_data['Times'])
x5 = toarray(cloud_frac_data['Times'])

y1 = toarray(Aerosol_data['Anomalies'])
y2 = toarray(NDVI_data['Anomalies'])
y3 = toarray(Precipitation_data['Anomalies'])
y4a = toarray(PAR_Total_data['Anomalies'])
y4b = toarray(PAR_Diffuse_data['Anomalies'])
y5 = toarray(cloud_frac_data['Anomalies'])

count = 2002
#print count + 0.875, count + 0.625



p = gridplot(x1,x2,x3,x4a,x4b,x5,y1,y2,y3,y4a,y4b,y5,figpath,"monthlyanomalies_r.pdf")
print p

#--------------------------------------------------
#AOD anomalie vs other anomalies
    
AOD_m = Aerosol_data['AOD Monthly Average']
Precip_m = Precipitation_data['Average Precipitation']
NDVI_m = NDVI_data['Average NDVI']
cloud_frac_m = cloud_frac_data['cloud_frac Monthly Average']
PAR_Total_m = PAR_Total_data['PAR Monthly Average']
PAR_Diffuse_m = PAR_Diffuse_data['PAR Monthly Average']

#plt.plot(Aerosol_data['AOD Monthly Average'], Precipitation_data['Average Precipitation'],'o')

#rainy season in the amazon occurs between December and May
#Dry season is from June to November




def season_sorting(x,y):
    wet = []
    dry = []
    for i in range(len(x)): 
        if (x[i] %1.0) > 0.875 or (x[i] %1.0) < 0.625: # rainy season months
            wet.append(i)
        if (x[i] %1.0) >= 0.625 and (x[i] %1.0) <= 0.875: # dry season months 
            dry.append(i)
    ywet = y[wet]
    ydry = y[dry]
    xwet = x[wet]
    xdry = x[dry]
    out = {'rain season x' : xwet, 'dry season x': xdry, 'rain season y': ywet, 'dry season y': ydry}
    return out

#aerosol
y1wet = season_sorting(x1,y1)['rain season y']
x1wet = season_sorting(x1,y1)['rain season x']

#ndvi
y2wet = season_sorting(x2,y2)['rain season y']
x2wet = season_sorting(x2,y2)['rain season x']

#precipitation
y3wet = season_sorting(x3,y3)['rain season y']
x3wet = season_sorting(x3,y3)['rain season x']
y3dry = season_sorting(x3,y3)['dry season y']
x3dry = season_sorting(x3,y3)['dry season x']

#par_total
#y4awet = 
#y4awey = 

#cloud fraction
y5wet = season_sorting(x5,y5)['rain season y']
x5wet = season_sorting(x5,y5)['rain season x']




cf_at_aerosol_wet = y5wet[xtime_equals_ytime(x5wet, x1wet)]
ndvi_at_cf_wet = y2wet[xtime_equals_ytime(x2wet, x5wet)]
ndvi_at_aero_wet = y2wet[xtime_equals_ytime(x1wet, x2wet)]
precip_at_cf_wet = y3wet[xtime_equals_ytime(x3wet,x5wet)]

print precip_at_cf_wet

lineofbestfit(ndvi_at_aero_wet,y1wet, 'Aero_rainy', 'NDVI_rainy',figpath, 'aero_ndvi_rainy.pdf')
lineofbestfit(y5wet, ndvi_at_cf_wet, 'cf wet', 'ndvi wet', figpath, 'cf_ndvi_rainy.pdf')
lineofbestfit( cf_at_aerosol_wet,y1wet, 'aerosol wet', 'cf wet', figpath, 'cf_aerosol_rainy.pdf')
lineofbestfit(precip_at_cf_wet,y5wet, ' precip wet', 'cloud_wet', figpath, 'precip_cloud_rainy.pdf')




sys.exit()

print precips, x3,x5
#precipitation vs cloud fraction
precip_at_cf = precips[xtime_equals_ytime(x3,x5)]
lineofbestfit(precip_at_cf, cfs, 'precip', 'cloudfrac',figpath, 'cfvprecip.pdf')

#NDVI vs PAR
ndvi_at_PAR = ndvis[xtime_equals_ytime(x2,x4b)]
lineofbestfit(ndvi_at_PAR,PAR_d, 'NDVI','PAR_diffuse',figpath,'ndvi_PAR_d.pdf')

#NDVI vs CF
ndvi_at_cf = ndvis[xtime_equals_ytime(x2,x5)]
lineofbestfit(ndvi_at_cf,cfs,'NDVI','Cloud Fraction', figpath, 'ndvi_cd.pdf')

#AOD vs Precipitation
precip_at_AOD = precips[xtime_equals_ytime(x1,x3)]
lineofbestfit(precip_at_AOD,aods, 'Precipitation', 'AOD', figpath, 'precip_aod.pdf')
sys.exit()
'''
y2aty1 = y2[xtime_equals_ytime(x1,x2)] # ndvi
y3aty1 = y3[xtime_equals_ytime(x1,x3)] # precipitation
y4aaty1 = y4a[xtime_equals_ytime(x1,x4a)] # PAR_total
y4baty1 = y4b[xtime_equals_ytime(x1,x4b)] # PAR_diffuse 
y5aty1 = y5[xtime_equals_ytime(x1,x5)] # Cloud Fraction
y2aty5 = y2[xtime_equals_ytime(x5,x2)] # NDVI when cloud frac readings were taken
y2aty4a = y2[xtime_equals_ytime(x4a,x2)] # NDVI when Total PAR readings were taken
y2aty4b = y2[xtime_equals_ytime(x4b,x2)] # NDVI when Diffuse PAR readings were taken 
'''


'''
y3aty5 = y3[xtime_equals_ytime(x5,x3)]

lineofbestfit(y2aty1,y1,'NDVI', 'Aerosol', 'NDVI_Aero_r.pdf')
lineofbestfit(y3aty1,y1,'Precipitation','Aerosol', 'Precip_Aero_r.pdf')
lineofbestfit(y5aty1,y1,'Cloud Fraction','Aerosol', 'Cloud_Aero_r.pdf' )
lineofbestfit(y3aty5,y5,'Precipitation Rate (mm/hr)', 'Cloud Fraction', 'Precip_Cloud_r.pdf')
lineofbestfit(x1,y1, 'Time','NDVI','NDVI_Time_r.pdf' )
lineofbestfit(y2aty5,y5,'NDVI','Cloud Fraction','NDVI_Cloud_r.pdf')
lineofbestfit(y2aty4a,y4a, 'NDVI','PAR Total (W/m2)', 'NDVI_PAR_Tot_r.pdf')
lineofbestfit(y2aty4b,y4b, 'NDVI','PAR Diffuse (W/m2)', 'NDVI_PAR_Dif_r.pdf')
'''





Average_monthlyNDVI = scatter(months,xticks,NDVI_data['Average NDVI'],'Average NDVI', 'Time(Monthly)', 'Average Monthly NDVI 2000-2017 - Rondonia (10S,64W) \n','monthly_NDVI_mean_plot_r.pdf', figpath,'xkcd:green') 
Average_monthly_cloud_frac = scatter(months,xticks,cloud_frac_data['cloud_frac Monthly Average'], 'Average Cloud Fraction','Time(Monthly)','Average Cloud Fraction per Month, July 2000 - April 2012 - Rondonia (10S,64W) \n','cloud_frac_monthly_average_r.pdf',figpath, 'xkcd:deep pink')
Average_monthlyPrecip = scatter(months,xticks,Precipitation_data['Average Precipitation'], 'Average Precipitation (mm/hr)','Time(Monthly)','Average Precipitation per Month, May 2000- March 2017 - Rondonia (10S,64W) \n','precip_monthly_average_r.pdf', figpath,'xkcd:violet')
Average_monthlyAOD = scatter(months,xticks,Aerosol_data['AOD Monthly Average'],'Average AOD', 'Time(Monthly)', 'Average Monthly AOD 2002-2012 - Rondonia (10S,64W) \n','monthly_AOD_mean_plot_r.pdf', figpath,'xkcd:sky blue')
Average_monthly_PAR_Total = scatter(months,xticks,PAR_Diffuse_data['PAR Monthly Average'], 'Average Flux (W/m2)','Time(Monthly)','Average Base of Atmosphere PAR (Diffuse) Flux per \n Month, July 2002 - April 2012 - Rondonia (10S,64W) \n','PAR_diffuse_monthly_average_r.pdf',figpath, 'xkcd:brown')
Average_monthly_PAR_Diffuse = scatter(months,xticks,PAR_Total_data['PAR Monthly Average'], 'Average Flux (W/m2)','Time(Monthly)','Average Base of Atmosphere PAR (Total) Flux per \n Month, July 2002 - April 2012 - Rondonia (10S,64W) \n','PAR_total_monthly_average_r.pdf',figpath, 'xkcd:tangerine')


'''
MonthlyAOD = scatter_reshape_axis(Aerosol_data['Times'],Aerosol_data['AOD Measurements'],year_list_AOD,'AOD', 'Time(Monthly)','Monthly AOD 2002-2012 - Rondonia (10S,64W) \n','monthly_AOD_2002_2012_r.pdf')
#AOD_num_retrievals = scatter(months,xticks,Aerosol_data['Monthly Retrievals'],'Monthly Retrievals','Time(Monthly)','Total Number of Retrievals - Rondonia (10S,64W) \n','AOD_monthly_retrievals_plot_r.pdf', 'xkcd:sky blue')
#Average_monthlyAOD = scatter(months,xticks,Aerosol_data['AOD Monthly Average'],'Average AOD', 'Time(Monthly)', 'Average Monthly AOD 2002-2012 - Rondonia (10S,64W) \n','monthly_AOD_mean_plot_r.pdf', 'xkcd:sky blue')
AOD_Anomalies_LOBF = scatter_LOBF(Aerosol_data['Times'],Aerosol_data['Anomalies'], year_list_AOD,'AOD Anomaly', 'Time(Monthly)', 'AOD Anomalies', 'AOD anomalies 2002-201 - Rondonia (10S,64W) \n','AOD_anomalies_rondonia_r.pdf'  )

MonthlyPrecip = scatter_reshape_axis(Precipitation_data['Times'], Precipitation_data['Precipitation'],year_list_PRECIP,'Mean Precipitation(mm/hr)','Time(Monthly)', 'Average Precipitation, May 2000- March 2017 \n - Rondonia (10S,64W) \n','monthly_AOD_2002_2012_r.pdf')
Precip_num_retreivals = scatter(months,xticks,Precipitation_data['Monthly Retrievals'],'Total Number of Retrievals','Time(Monthly)','Monthly Retrievals - Precipitation - Rondonia (10S,64W)\n','monthly_retrievlas_precip_r.pdf', 'xkcd:voilet')
#Average_monthlyPrecip = scatter(months,xticks,Precipitation_data['Average Precipitation'], 'Average Precipitation (mm/hr)','Time(Monthly)','Average Precipitation per Month, May 2000- March 2017 - Rondonia (10S,64W) \n','precip_monthly_average_r.pdf')
Precip_Anomalies_LOBF = scatter_LOBF(Precipitation_data['Times'],Precipitation_data['Anomalies'],year_list_PRECIP, 'Precipitation Anomaly (mm/hr)','Time(Yearly)','Precipitation Anomalies','Precipitation Anomalies May 2000 - March 2017 - Rondonia (10S,64W) \n','Precipitation_anomalies_r.pdf' )

#MonthlyPAR = scatter_reshape_axis(PAR_data['Times'], PAR_data['PARs'],year_list_AOD,'Mean top of atmosphere \n shortwave downwelling flux (W/m^2)','Time(Monthly)', 'Average shortwave downwelling flux,\n July 2000 - April 2012 - Rondonia (10S,64W) \n','monthly_PAR_2002_2012_r.pdf')
#PAR_Diffuse_num_retreivals = scatter(months,xticks,PAR_data['Monthly Retrievals'],'Total Number of Retrievals','Time(Monthly)','Monthly Retrievals - Top of atmosphere \n shortwave downwelling flux - Rondonia (10S,64W) \n','monthly_retrievlas_PAR_r.pdf', 'xkcd:sky blue')
#Average_monthly_PAR = scatter(months,xticks,PAR_data['PAR Monthly Average'], 'Average top of atmosphere \n shortwave downwelling flux (W/m^2)','Time(Monthly)','Average top of atmosphere shortwave downwelling \n flux per Month, July 2000 - April 2012 - Rondonia (10S,64W) \n','PAR_monthly_average_r.pdf')
#PAR_Anomalies_LOBF = scatter_LOBF(PAR_data['Times'],PAR_data['Anomalies'],year_list_AOD, 'Top of atmosphere shortwave \n downwelling flux anomaly (W/m^2)','Time (Yearly)','Top of atmosphere shortwave \n downwelling flux anomalies','Top of atmosphere shortwave downwelling flux \n aomalies, July 2000 - April 2012 - Rondonia (10S,64W) \n','PAR_anomalies_r.pdf' )

Monthlycloud_frac = scatter_reshape_axis(cloud_frac_data['Times'], cloud_frac_data['cloud_fracs'],year_list_AOD,'Mean cloud fraction','Time(Monthly)', 'Average cloud fraction \n July 2000 - April 2012 - Rondonia (10S,64W) \n','monthly_cloud_frac_2002_2012_r.pdf')
#cloud_frac_num_retreivals = scatter(months,xticks,cloud_frac_data['Monthly Retrievals'],'Total Number of Retrievals','Time(Monthly)','Monthly Retrievals - cloud fraction - Rondonia (10S,64W) \n','monthly_retrievlas_cloud_frac_r.pdf')
#Average_monthly_cloud_frac = scatter(months,xticks,cloud_frac_data['cloud_frac Monthly Average'], 'Average cloud fraction','Time(Monthly)','Average cloud fraction per Month, July 2000 - April 2012 - Rondonia (10S,64W) \n','cloud_frac_monthly_average_r.pdf')
cloud_frac_Anomalies_LOBF = scatter_LOBF(cloud_frac_data['Times'],cloud_frac_data['Anomalies'],year_list_AOD, 'cloud fraction anomaly','Time (Yearly)','Cloud fraction anomalies','cloud fraction aomalies, July 2000 - April 2012 - Rondonia (10S,64W) \n','cloud_frac_anomalies_r.pdf' )

MonthylyNDVI = scatter_reshape_axis(NDVI_data['Times'], NDVI_data['NDVIs'], year_list_PRECIP, 'NDVIs', 'Time(Monthly)', 'Monthly NDVI 2000-2017 - Rondonia (10S,64W) \n','monthly_NDVI_2002_2012_r.pdf')
#NDVI_num_retrievals = scatter(months,xticks,NDVI_data['Monthly Retrievals'],'Monthly Retrievals','Time(Monthly)','Total Number of NDVI Retrievals - Rondonia (10S,64W) \n','NDVI_monthly_retrievals_plot_r.pdf')
#Average_monthlyNDVI = scatter(months,xticks,NDVI_data['Average NDVI'],'Average NDVI', 'Time(Monthly)', 'Average Monthly NDVI 2000-2017 - Rondonia (10S,64W) \n','monthly_NDVI_mean_plot_r.pdf')
NDVI_Anomalies_LOBF = scatter_LOBF(NDVI_data['Times'],NDVI_data['Anomalies'], year_list_PRECIP,'NDVI Anomaly', 'Time(Monthly)', 'NDVI Anomalies', 'NDVI anomalies 2000-2017 - Rondonia (10S,64W) \n','NDVI_anomalies_rondonia.pdf'  )
'''