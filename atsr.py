
#---------------------------------------------------------------------------------------------
# Name: Module containing ATSR file read functions
# Functions: aerosol_file_info, aerosol_file_attributes
# History:
#   07/27/17 MC: add ATSR filename reader module for aerosol files
#---------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------
# Name: aerosol_file_info
# Input: fname: name of the AEROSOL file with path included
# Output: tOUT: dictionary containing the filename without path, month and year
#---------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

def aerosol_file_info(fname):
    # seperating the file path by '/'
    getfilepath = fname.split('/')
    # Date and month info is found in the last item in this list
    filename = getfilepath[(len(getfilepath)-1)]
    #Extract the year
    YYYY = filename[:4]
    YEAR = float(YYYY)
    #Extract the month
    MM = filename[4:6]
    #print MM
    MONTH = float(MM)
    OUT = { "filename" : filename, "MONTH" : MONTH, "MM" :MM, "YEAR" : YEAR, "YYYY" : YYYY }
    return OUT

#---------------------------------------------------------------------------------------------
# Name: fetch_aerosol_file_attributes
# Description: essentially runs the above function but for a list of files (not just one)
# Input: list_of_files (files from which time will be extracted)
# Output: dictionary containing all of the files and times
# Improvements: add cloud reader, could output file type, julian day, ect...
#---------------------------------------------------------------------------------------------
def fetch_aerosol_file_attributes(list_of_files):
 
    list_of_times = []
    list_of_months = [] 
    list_of_years = []
    for item in list_of_files:
        item_attributes = aerosol_file_info(item)
        list_of_months.append(item_attributes['MONTH'])
        list_of_years.append(item_attributes['YEAR'])
        ts = item_attributes['YEAR'] + (item_attributes['MONTH']-0.5)/12. #centre of the month
        list_of_times.append(ts)# writing each time to a new list
    tOUT= { "time" : list_of_times, "month" : list_of_months, "year" : list_of_years, "file" : list_of_files }
    return tOUT


#--------------------------------------------------------------------------------------------

def landuse_file_info(fname):

    # seperating the file path by '/'
    getfilepath = fname.split('/')

    # Date and month info is found in the last item in this list
    filename = getfilepath[(len(getfilepath)-1)] # last element in list
    #Extract the year
    YYYY = filename[:4]
    YEAR = float(YYYY)
    #Extract the month
    MM = filename[4:6]
    #print MM
    MONTH = float(MM)

    OUT = { "filename" : filename, "MONTH" : MONTH, "MM" :MM, "YEAR" : YEAR, "YYYY" : YYYY }
    return OUT

#------------------------------------------------------------------------------------------------


def fetch_landuse_file_attributes(list_of_files):
 
    list_of_times = []
    list_of_months = [] 
    list_of_years = []
    for item in list_of_files:
        item_attributes = landuse_file_info(item)
        list_of_months.append(item_attributes['MONTH'])
        list_of_years.append(item_attributes['YEAR'])
        ts = item_attributes['YEAR'] + (item_attributes['MONTH']-0.5)/12. #centre of the month
        list_of_times.append(ts)# writing each time to a new list
    tOUT= { "time" : list_of_times, "month" : list_of_months, "year" : list_of_years, "file" : list_of_files }
    return tOUT

    
def seasonal_anomalies(Years, Months, Times, Data, datatitledry,datatitlewet, dataname):
    ti = 2003
    rainy_season_averages = []
    dry_season_averages = []
    rainyID = []
    dryID = []
    rainy_times = []
    dry_times = []
    wet_anoms = []
    dry_anoms = []
    num_retrivals_wet = []
    num_retrivals_dry = []
    
    for i in range(9):
        YID_rain = np.where((Years == ti) & ((Months == 12) | (Months < 6)))
        rainyID.append(YID_rain)
        wet_anoms.extend(Data[YID_rain])
        num_retrivals_wet.append((Data[YID_rain]).size)
        yearly_mean_rain = np.mean(Data[YID_rain]) # calculating the mean for each year for the rainy season
        rainy_times.extend(Times[YID_rain])
        rainy_season_averages.append(yearly_mean_rain)
        YID_dry = np.where((Years == ti) & ((Months > 5) & (Months < 12))) 
        dryID.append(YID_dry)
        dry_anoms.extend(Data[YID_rain])
        dry_times.extend(Times[YID_dry])
        yearly_mean_dry = np.mean(Data[YID_dry]) # calculating the mean for each year for the dry season
        dry_season_averages.append(yearly_mean_dry)
        num_retrivals_dry.append((Data[YID_dry]).size)
        ti = ti+1 
        
    nine_yr_mean_r = np.mean(wet_anoms)
    nine_yr_mean_d = np.mean(dry_anoms)
    
    
    nine_yr_anom_wet = []
    for i in range(9):
        r_anom = rainy_season_averages[i] - nine_yr_mean_r
        nine_yr_anom_wet.append(r_anom)
    
    nine_yr_anom_dry = []
    for i in range(9):
        d_anom = dry_season_averages[i] - nine_yr_mean_d
        nine_yr_anom_dry.append(d_anom)
    # calculating anomalies by subtrating the seasonal means 
    r_anoms = []
    d_anoms = []    
    for i in range(len(rainyID)): 
        r_anomaly = Data[rainyID[i]] - rainy_season_averages[i]
        r_anoms.extend(r_anomaly)
        
    for i in range(len(dryID)):
        d_anomaly = Data[dryID[i]] - dry_season_averages[i]
        d_anoms.extend(d_anomaly)
    '''     
    plt.plot(np.linspace(2003,2011,9), num_retrivals_dry)
    plt.title(datatitledry)
    plt.show()
    
    plt.plot(np.linspace(2003,2011,9), num_retrivals_wet)
    plt.title(datatitlewet)
    plt.show()
    
    plt.plot((np.linspace(2003,2011,9)), dry_season_averages)
    plt.title(datatitledry)
    plt.xlabel('time (year)')
    plt.ylabel(dataname)
    plt.show()
    '''
    out = {'dry anomalies': d_anoms, 'rainy anomalies': r_anoms, 'dry times': dry_times, 'rainy times': rainy_times, 'nine year anomalies rain': nine_yr_anom_wet, 'nine year anomalies dry': nine_yr_anom_dry}
    return out