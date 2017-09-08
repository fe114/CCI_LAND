#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 13:26:03 2017

@author: fespir
"""

def PAR_file_info(fname):
    # seperating the file path by '/'
    getfilepath = fname.split('/')
    # Date and month info is found in the last item in this list
    filename = getfilepath[(len(getfilepath)-1)]
    #print filename
    namelist = filename.split('_')
    yr_m = namelist[2]# last element in list
    #Extract the year
    YYYY = yr_m[:4]
    YEAR = float(YYYY)
    #Extract the month
    MM = yr_m[4:6]
    #print MM
    MONTH = float(MM)
    OUT = { "filename" : filename, "MONTH" : MONTH, "YEAR" : YEAR}
    return OUT

def get_PAR_attributes(filelist):
    month = []
    year = []
    time = []
    files = []
    for item in filelist:
        extractdata = PAR_file_info(item)
        M = float(extractdata["MONTH"])
        Y = float(extractdata["YEAR"])
        month.append(M)
        year.append(Y)
        time.append(Y+((M-0.5)/12))
        files.append(item)
    dictionary = {"Year": year, "Month": month, "Time": time, "File": files}
    return dictionary
        
      