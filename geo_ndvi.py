#Procedure to obtain latitude and longitude values for
#MODIS NDVI data

from numpy import *
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

    OUT = { "latitude" : lat, "longitude" : lon }
    return OUT
