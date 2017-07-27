#---------------------------------------------------------------------------------------------
# Name: Module containing ATSR functions
#
# Input: fname: name of the file with path included for standard aerosol files
#
# Output: tOUT: dictionary containing the filename without path, month and year
#
# Improvements: add cloud reader, could output file type, julian day, ect...
#
# History:
#   07/27/17 MC: add ATSR filename reader module for aerosol files
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

    tOUT = { "filename" : filename, "MONTH" : MONTH, "MM" :MM, "YEAR" : YEAR, "YYYY" : YYYY }
    return tOUT
