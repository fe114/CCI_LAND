#---------------------------------------------------------------------------------------------
# Name: file_search module
#
# Input: top-level path to files, suffix files ending with (e.g. ".nc")        
#
# Output: all of the files
#
# History:
#   07/31/17 MC: add file_search.py module
#---------------------------------------------------------------------------------------------
import os,sys
def extract_files(path,suffix):

    list_of_paths = []
    for root, dirs, files in os.walk(path): 
        for file in files:
            if file.endswith(suffix):
                file_path = (os.path.join(root, file))
                list_of_paths.append(file_path) # writes each full file path to the empty list defined above
    
    return list_of_paths
