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
def aerosol_file_info(fname):

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

