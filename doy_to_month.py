#-------------------------------------
# Name: convert_doy_to_month
# Input: doy
# Output: month
# module: jdcal (e.g. conda install jdcal)
#-------------------------------------
from jdcal import *
def doy_to_month(doy,year):
    #fetch julian day
    jday0 = gcal2jd(year,1,1)
    jday1 = (jday0[0],jday0[1]+doy-1)
    caldat = jd2gcal(*jday1)
    month = caldat[1]
    return month
