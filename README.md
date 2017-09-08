# CCI_LAND
# Python code to process Aerosol (.nc), Precipitation(.nc), Cloud Fraction (.nc), PAR(.nc) and NDVI (.hdf) files
# Framework is set up for jasmin-cems 
 

User guide: 
Process data for a lat:lon coordinate grid of interest. The grid must be of equal dimensions e.g. for Rondonia, 
lat boundary = [-11,-9], lon boundary = [-65,-63]. Processing is done in top_level_code.py
1. Open top_level_code.py and change the outpath and figpaths 
2. Run the program and input the lat and lon upper and lower boundaries when prompted
3. You will be asked if you want to remove any previously made files which were made if the program was run before using the same
lat lon region. Print y when prompted to delete them or n to keep them. 
This is necessary if you want to reprocess data for the same region as the savefiles will not be reproduced if the lat lon boundaries are constant 
4. Once processed (takes ~ 30 mins), open data2plots.py
5. Change the outpath and figpath to the same paths used in top_level_code.py
6. Run the program, enter the lat lon boundary for the files produces by top_level_code.py and a variety of plots will appear! 
