#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 10:39:17 2017

@author: fespir
"""


from numpy import * 
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
import csv

figpath = '/group_workspaces/cems2/nceo_generic/CCI_LAND/figs/'
aerosol_path = "/group_workspaces/cems/aerosol_cci/public/cci_products/AATSR_ORAC_v04-01/L3_MONTHLY/"
precip_path = "/group_workspaces/cems2/nceo_generic/satellite_data/precipitation/TRMM/monthly"
ndvi_path = "/group_workspaces/cems2/nceo_generic/satellite_data/modis_c6/ndvi/"
outpath = "/group_workspaces/cems2/nceo_generic/CCI_LAND/output_data/"

lon = -63.5 
lat = -11.5 
#lat_bnds, lon_bnds = [-12, -11], [-64, -63]
lat_bnds, lon_bnds = [-11, -9], [-65.5, -63.5]

Aerosol_data = process_aerosol(lon,lat,outpath,aerosol_path)
Precipitation_data = process_precip(lon_bnds,lat_bnds,outpath,precip_path) 
NDVI_data = process_ndvi(lon_bnds,lat_bnds,outpath,ndvi_path,".hdf")
#Landuse_data = 


print NDVI_data

sys.exit()







#--------------------------------
#plotting

xticks = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'] 
year_list_AOD = linspace(2002,2012,11)
months = linspace(1,12,12)
year_list_PRECIP = linspace(2000,2018,19)


def scatter_reshape_axis(Time, Data,xticks, ylab, xlab, figtitle, figname):
    #plot AOD readings for each month
    fig = plt.figure(figsize=(16,4))
    ax = fig.add_subplot(111)
    plt.plot(Time, Data, 'o-')
    plt.title(figtitle, fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xticks(xticks,fontsize = 14)
    plt.ylabel(ylab, fontsize = 14)
    plt.xlabel(xlab, fontsize = 14)
    figshow = plt.show()
    plt.savefig( figpath + figname)
    return figshow


def scatter(tick_spacing,xticks, retrievals, ylab, xlab, figtitle, figname):
    #number of retrievlas for each month from 2002-2012
    plt.xticks(tick_spacing,xticks, fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.plot(tick_spacing,retrievals, 'o-')
    plt.title(figtitle, fontsize = 14)
    plt.ylabel(ylab, fontsize = 14)
    plt.xlabel(xlab, fontsize = 14)
    plot = plt.show()
    plt.savefig(figpath+figname)
    return plot



def scatter_LOBF(Time, Data, xticks, ylab, xlab, scatterlabel, figtitle, figname):
#plotting AOD anomalies
    fig = plt.figure(figsize=(16,4))
    ax = fig.add_subplot(111)
    plt.title(figtitle, fontsize = 14)
    plt.xlabel(xlab,fontsize = 14 )
    plt.ylabel(ylab, fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xticks(xticks,fontsize = 14)
    plot = plt.plot(Time, Data, 'o', label=scatterlabel, markersize=5) 
    matrix = np.vstack([Time, np.ones(len(Time))]).T
    m, c = np.linalg.lstsq(matrix, Data)[0]

    #finding intercept at the time of the first AOD reading
    y_x0= m*Time[0] + c

    #.5f =  5 significant figures
    bestfittxt =  '{:.5f}x + {:.10f}'.format(m, y_x0)
    
    #plotting line of best fit 
    bestfit = plt.plot(Time, m*Time + c, 'r', label='Fitted line ' + bestfittxt)

    #calculating the mean, variance, T statistic and p value
    m, v, s, k = stats.t.stats(5, moments='mvsk')
    n, (smin, smax), sm, sv, ss, sk = stats.describe(Data)
    sstr = 'mean = %6.4f, variance = %6.4f'
    tt = (sm-m)/np.sqrt(sv/float(n))  # t-statistic for mean
    pval = stats.t.sf(np.abs(tt), n-1)*2  # two-sided pvalue = Prob(abs(t)>tt)
    print 't-statistic = %6.3f pvalue = %6.4f' % (tt, pval)
    print 'distribution: ', sstr %(m, v)
    print 'sample: ', sstr %(sm, sv)
    plt.legend()
    plt.show()
    plt.savefig( figpath + figname)
    return plot, bestfit 


MonthlyAOD = scatter_reshape_axis(Aerosol_data['Times'],Aerosol_data['AOD Measurements'],year_list_AOD,'AOD', 'Month','Monthly AOD 2002-2012 \n','monthly_AOD_2002_2012.pdf')
AOD_num_retrievals = scatter(months,xticks,Aerosol_data['Monthly Retrievals'],'Month','Monthly Retrievals','Total Number of Retrievals \n','AOD_monthly_retrievals_plot.png')
Average_monthlyAOD = scatter(months,xticks,Aerosol_data['AOD Monthly Average'],'AOD', 'Month', 'Average Monthly AOD 2002-2012 \n','monthly_AOD_mean_plot.png')
AOD_Anomalies_LOBF = scatter_LOBF(Aerosol_data['Times'],Aerosol_data['Anomalies'], year_list_AOD,'AOD Anomaly', 'Time (month)', 'AOD Anomalies', 'AOD anomalies 2002-2012, Rondonia \n','AOD_anomalies_rondonia.pdf'  )

MonthlyPrecip = scatter_reshape_axis(Precipitation_data['Times'], Precipitation_data['Precipitation'],year_list_PRECIP,'Mean Precipitation(mm/hr)','Month', 'Average Precipitation, May 2000- March 2017 \n','monthly_AOD_2002_2012.pdf')
Precip_num_retreivals = scatter(months,xticks,Precipitation_data['Monthly Retrievals'],'Total Number of Retrievals','Month','Monthly Retrievals - Precipitation \n','monthly_retrievlas_precip.png')
Average_monthlyPrecip = scatter(months,xticks,Precipitation_data['Average Precipitation'], 'Average Monthly Precipitation (mm/hr)','Month','Average Precipitation per Month, May 2000- March 2017 \n','precip_monthly_average.png')
recip_Anomalies_LOBF = scatter_LOBF(Precipitation_data['Times'],Precipitation_data['Anomalies'],year_list_PRECIP, 'Precipitation Anomaly (mm/hr)','Year','Precipitation Anomalies','Precipitation Anomalies May 2000 - March 2017 \n','Precipitation_anomalies.png' )