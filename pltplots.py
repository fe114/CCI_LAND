#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 16:24:23 2017

@author: fespir

NAME;
Program to produce plots and statistics 

PURPOSE;
statistics and list processing:
1. Moving average (movingaverage)
2. Line of best fit (LOBF) --> used for a gridplot e.g. 3x2
3. Best Fit Line (bestfitline) very similar to 2. but used for one single plot
4. t statistic (tstat)
5. Converts two different length arrays into two equal length arrays (xtime_equals_ytime)
6. list to array converter (toarray)


Plots:
1. 3x1 plot, x axis is shared by the three plots (triplotsharex)
2. 3x2 grid plot, one column of three plots shares one x axis (gridplot). Moving average funtion is applied to each plot. 
3. scatter plot using two lists (scatter)
4. scatter plot with line of best fit (lineofbestfit). T statsic is ploted in legend
5. scatter plot which is not saved to a file  


AUTHOR;
Freya Espir

"""

from numpy import * 
import math
import os,sys
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt 
from scipy import stats 
import csv
from matplotlib import gridspec
from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go
from scipy.stats import linregress

#-----------------------------------


def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

def LOBF(x,y,axes1,axes2,ax):
    matrix = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(matrix, y)[0]
    #finding intercept at the time of the first AOD reading
    y_x1= m*x[0] + c
    bestfittxt =  '{:.3f}x + {:.3f}'.format(m, y_x1)
    #plotting line of best fit 
    bestfit = ax[axes1,axes2].plot(x, m*x+ c, 'r',label= 'LOBF ' + bestfittxt)
    ax[axes1,axes2].legend()
    return bestfit

def bestfitline(x,y):
    matrix = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(matrix, y)[0]

    #finding intercept at the time of the first AOD reading
    x_1 = m*x[0] + c

    bestfittxt =  '{:.3f}x + {:.3f}'.format(m, x_1)
    slope = '{:.3f}'.format(m)
    #plotting line of best fit 
    bestfit = plt.plot(x, m*x+ c, label=bestfittxt)
    r = linregress(x,y)[2]
    N = x.size
    Ttest = (r*np.sqrt(N))/(np.sqrt(1.0-(r**2)))
    #tstat ='t-stat =  ' + '{:.3f}'.format(Ttest)# t-statistic for mean
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    return bestfittxt

def bestfitlinemultiplot(x,y,axnum): 
    matrix = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(matrix, y)[0]
    #finding intercept at the time of the first AOD reading
    y_x1= m*x[0] + c
    bestfittxt =  '{:.2f}x + {:.2f}'.format(m, y_x1)
    #plotting line of best fit 
    bestfit = axnum.plot(x, m*x+ c, 'r',label= bestfittxt)
    
    axnum.legend()
    return bestfit
    
def tstat(x,y):
    r = linregress(x,y)[2]
    N = x.size
    Ttest = (r*np.sqrt(N))/(np.sqrt(1.0-(r**2)))
    tstat ='t-stat =  ' + '{:.3f}'.format(Ttest)# t-statistic for mean
    return tstat

def triplotsharex(x1,x2,x3,y1,y2,y3,xlab,ax1name,ax2name,ax3name,path,figname,figtitle):
    f,(ax1,ax2,ax3) = plt.subplots(3,1,sharex= True)
    ax1.plot(x1,y1,'o',color = 'xkcd:sky blue', label = ax1name)
    plot1bestfit = bestfitlinemultiplot(x1,y1,ax1)
    
    ax2.plot(x2, y2, 'o',color ='xkcd:green', label = ax2name)
    plot2bestfit = bestfitlinemultiplot(x2,y2,ax2)
    
    ax3.plot(x3,y3,'o',color = 'xkcd:violet',label = ax3name)
    plot3bestfit = bestfitlinemultiplot(x3,y3,ax3)
    
    #changing axis size whilst keeping scale constant
    F = plt.gcf()
    DefaultSize = F.get_size_inches()
    F.set_size_inches( (DefaultSize[0], DefaultSize[1]*2) )
    Size = F.get_size_inches()
    plt.subplots_adjust(hspace=.1)
    plt.xlabel(xlab)
    f.text(0.3,0.9, figtitle, va='center', rotation='horizontal', fontsize = 12 )
    plt.show()
    f.savefig(path + figname,bbox_inches='tight')
   
    
def gridplot(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,path,figname,figtitle):
    
    f, axarr = plt.subplots(3,2, sharex=True)
    Aerosol, = axarr[0,0].plot(x1, y1, 'xkcd:sky blue',label='Aerosol')
    y1_av = movingaverage(y1, 12)
    y1ave, = axarr[0,0].plot(x1, y1_av, "--", color = "k")#, label = 'yearly running mean')
    y1bestfit = LOBF(x1,y1,0,0,axarr)
    f.text(0.12,0.89, tstat(x1,y1), color = 'r', va='center', rotation='horizontal', fontsize = 10 )
    
    NDVI, = axarr[0,1].plot(x2, y2, color ='xkcd:green',label='NDVI')
    y2_ave = movingaverage(y2, 12)
    y2ave, = axarr[0,1].plot(x2,y2_ave, "--",color = "k")#label = 'yearly running mean')
    y2bestfit = LOBF(x2,y2,0,1,axarr)
    f.text(0.56,0.89, tstat(x2,y2), color = 'r', va='center', rotation='horizontal', fontsize = 10 )
    
    Precip, = axarr[1,0].plot(x3, y3, color = 'xkcd:violet', label = 'Precipitation Rate (mm/hr)')
    y3_ave = movingaverage(y3,12)
    y3ave, = axarr[1,0].plot(x3,y3_ave, "--",color = "k")#, label = 'yearly running mean')
    y3bestfit = LOBF(x3,y3,1,0,axarr)
    f.text(0.12,0.63, tstat(x3,y3), color = 'r', va='center', rotation='horizontal', fontsize = 10 )
    
    PAR_total, = axarr[1,1].plot(x4, y4, color = 'xkcd:tangerine', label = 'PAR_total(W/m2)')
    y4_ave = movingaverage(y4,12)
    y4ave, = axarr[1,1].plot(x4,y4_ave, "--",color = "k")#, label = 'yearly running mean')
    y4bestfit = LOBF(x4,y4,1,1,axarr)
    f.text(0.56,0.63, tstat(x4,y4), color = 'r', va='center', rotation='horizontal', fontsize = 10 )
    
    PAR_diffuse, = axarr[2,0].plot(x5, y5, color = 'xkcd:brown', label = 'PAR_diffuse (W/m2)')
    y5_ave = movingaverage(y5,12)
    y5ave, = axarr[2,0].plot(x5,y5_ave, "--",color = "k")#, label = 'yearly running mean')
    y5bestfit = LOBF(x5,y5,2,0,axarr)
    f.text(0.12,0.37, tstat(x5,y5), color = 'r', va='center', rotation='horizontal', fontsize = 10 )
    
    CCI, = axarr[2,1].plot(x6, y6, color = 'xkcd:deep pink', label = 'Cloud Fraction')
    y6_ave = movingaverage(y6,12)
    y6ave, = axarr[2,1].plot(x6,y6_ave, "--",color = "k")#, label = 'yearly running mean')
    y6bestfit = LOBF(x6,y6,2,1,axarr)
    f.text(0.56,0.37, tstat(x6,y6), color = 'r', va='center', rotation='horizontal', fontsize = 10 )
    
    #plt.xlabel('Time(Monthly)',fontsize = 16 )
    plt.xticks(linspace(2000,2018,19))
    #f.text(0.03,0.5, 'Anomalies', va='center', rotation='vertical', fontsize = 12 )
    f.text(0.45,0.08, 'Time (monthly)', va='center', rotation='horizontal', fontsize = 12 )
    f.text(0.35,0.92, figtitle, va='center', rotation='horizontal', fontsize = 12 )
    
    F = plt.gcf()
    DefaultSize = F.get_size_inches()
    print "Default size in Inches", DefaultSize
    F.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]*3) )
    Size = F.get_size_inches()
    print "Size in Inches", Size
    plt.subplots_adjust(hspace=.1)
    plt.xlim([2002,2012])

    #print fig_size
    
    #axarr[4].legend((Aerosol, NDVI, Precip, PAR, CCI), ('Aerosol', 'NDVI', 'Precipitation', 'PAR', 'CCI'), loc='lower left')

    plt.show()
    f.savefig(path + figname,bbox_inches='tight')
    


def scatter(tick_spacing,xticks, retrievals, ylab, xlab, figtitle, figname, path, plotcolour):
    fig = plt.figure()
    #number of retrievlas for each month from 2002-2012
    plt.xticks(tick_spacing,xticks, fontsize = 10)
    plt.yticks(fontsize = 10)
    plt.plot(tick_spacing,retrievals, 'o-', color = plotcolour)
    plt.title(figtitle, fontsize = 12)
    plt.ylabel(ylab, fontsize = 12)
    plt.xlabel(xlab, fontsize = 12)
    plot = plt.show()
    fig.savefig(path+figname,bbox_inches='tight')
    return plot


def toarray(mylist):
    return np.asarray(mylist)


def lineofbestfit(x,y,xlab,ylab,path,figname):
    #fig = plt.figure()
    fig = plt.figure(figsize=(4,4))
    x = np.asarray(x)
    y = np.asarray(y)
    plt.plot(x, y, 'o')
    matrix = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(matrix, y)[0]

    #finding intercept at the time of the first AOD reading
    x_1 = m*x[0] + c

    #.5f =  5 significant figures
    bestfittxt =  '{:.3f}x + {:.3f}'.format(m, x_1)
    slope = '{:.3f}'.format(m)
    #plotting line of best fit 
    bestfit = plt.plot(x, m*x+ c, 'r', label= 'slope =  ' + slope)
    r = linregress(x,y)[2]
    N = x.size
    Ttest = (r*np.sqrt(N))/(np.sqrt(1.0-(r**2)))
    tstat ='t-stat =  ' + '{:.3f}'.format(Ttest)# t-statistic for mean
    plt.xlabel(xlab, fontsize = 10)
    plt.ylabel(ylab, fontsize = 10)
    plt.legend()
    #fig.legend((theplot,bestfit,tt), (ylab,bestfit,tt))
    fig.text(0.15,0.9, tstat, color = 'r', va='center', rotation='horizontal', fontsize = 10 )
    plt.show()
    fig.savefig(path + figname, bbox_inches='tight')#, pad_inches=0)
    return bestfittxt





#lineofbestfit(AOD_m,NDVI_m)
#lineofbestfit(CCI_m, NDVI_m, 'Cloud Fraction monthly average', 'NDVI monthly average', )


def xtime_equals_ytime(axis1, axis2):
    location = []
    for i in range(axis1.size):
        for j in range(axis2.size):
            #print x1[i], x2[i]
            if axis1[i] == axis2[j]:
                point = i 
                location.append(point)
    loc = np.asarray(location)
    return loc


def t_test(x,y):
    r = linregress(x,y)[2]
    N = len(x)
    Ttest = (r*np.sqrt(N))/(np.sqrt(1.0-(r**2)))
    return Ttest, N, r, np.sqrt(N)


def scatterplot(x,y,xlab,ylab,figtitle):
    fig = plt.figure(figsize=(4,4))
    plt.plot(x,y,'o')
    plt.ylabel(ylab)
    plt.xlabel(xlab)
    plt.title(figtitle)
    plt.show()
    return fig