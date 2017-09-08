#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 10:00:20 2017

@author: fespir
"""

from solar import *

#rondonia lat and lon bounds

Lat_bounds, Lon_bounds = [-11,-9], [-65,-63]
latitude = (Lat_bounds[0] + Lat_bounds[1])/2.0

print latitude 
sys.exit()
solar_radiation_daily(latitude,caldat)