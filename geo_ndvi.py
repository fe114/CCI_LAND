from numpy import *
from doy_to_month import *
from pyhdf.SD import SD, SDC
import numpy as np

#Procedure to obtain latitude and longitude values for
#MODIS NDVI data
def geo_ndvi():

    varDim = 7200
    lon = zeros(varDim)
    ct=0
    for i in lon:
        lon[ct] = ct*360/float(varDim) - 180. + 0.025
        ct=ct+1
    varDim = 3600
    lat = zeros(varDim)
    ct=0
    for i in lat:
        lat[ct] = ct*180/float(varDim) - 90. + 0.025
        ct=ct+1
    lat = lat[::-1]
    OUT = { "latitude" : lat, "longitude" : lon }
    return OUT

#procedure to extract the year and day of year 
#from an .hdf file name and convert to month
def NDVI_file_info(fname):
    getfilepath = fname.split('.')
    filedate = getfilepath[1]
    YEAR = filedate[1:5]
    DAY = filedate[5:8]
    day2month = doy_to_month(int(DAY),int(YEAR))
    OUT = {"YEAR": int(YEAR),"DAY":int(DAY), "MONTH":day2month}
    return OUT

#write the NDVI_file_info output to lists
#create a list of times, where TIME = YEAR + (MONTH-0.5)/12
def fetch_NDVI_file_attributes(filelist):
    daylist = []
    yearlist = []
    month = []
    times = []
    for item in filelist:
        filename_attributes = NDVI_file_info(item)
        yearlist.append(filename_attributes['YEAR'])
        daylist.append(filename_attributes['DAY'])
        month.append(filename_attributes['MONTH'])
        ts = filename_attributes['YEAR'] + ((filename_attributes['MONTH']-0.5)/12.0)
        times.append(float(ts))
    OUT = {"Year":yearlist, "Month": month, "Day":daylist, "Times":times}
    return OUT #yearlist,month,daylist,times,files

#get the NDVI readings for a coodinate grid from the .hdf file
def get_NDVI_readings(item,lat_bnds,lon_bnds,datasetname,scale_factor,add_offset):
    file = SD(item, SDC.READ)
    datasets_dic = file.datasets()
    sds_obj = file.select(datasetname) # select sds
    geodata = geo_ndvi()
    lat = geodata['latitude'] #list of latitudes 
    lon = geodata['longitude'] 
    # list of longitudes
    lat_inds = np.where((lat > lat_bnds[0]) & (lat < lat_bnds[1]))[0] #index value of all latitudes within coodinate grid
    lon_inds = np.where((lon > lon_bnds[0]) & (lon < lon_bnds[1]))[0] #index value of all longiudes within coodinate grid
    lat_zero = lat_inds[0] # lower boundary latitude value
    lat_max = lat_inds[len(lat_inds)-1] # upper boundary latitude value
    lon_zero = lon_inds[0]  # lower boundary longitude value
    lon_max = lon_inds[len(lon_inds)-1]  # lower boundary longitude value
    data = sds_obj[lat_zero:lat_max,lon_zero:lon_max]#getting ndvi readings within lat,lon boundaries
    data = (data - add_offset) / scale_factor # scaling the data using scale_factor from .hdf file
    return data