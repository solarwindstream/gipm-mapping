##OMNI 4MINS REDO

import pandas as pd
import numpy as np
import glob

from omni_seg_centered_4mins import omni_seg

#OMNI CSVs made using CDF to CSV script

#omni_seg(om_df, only_full_windows): takes list of window start times
#and produces OMNI averages for list of windows (return(om_averages))
#maybe would be best off though doing omni_seg & attaching results TO CLUSTER FILES
#keep omni csvs too though as they're useful
#do omni_seg first, save them, & then pair matching CSVs and add relevant OMNI stats

#load Cluster GIPM files for 2001 (test), saving just their interval starts, and then do only six of them

list_all = []

path = "/Users/apx059/Documents/4_Min_GIPM/2022_Integrated_CSVs/**"

for path in glob.glob(path, recursive=True):
    list_all.append(path)
    
#list with only files, not folders
cl_csvs = []

for element in list_all:
    if '.csv' in element:
        cl_csvs.append(element)
        
cl_interval_starts = []

for element in cl_csvs:
    cl = pd.read_csv(element)
    cl['datetime'] = pd.to_datetime(cl['datetime'],format='mixed')
    interval_list = cl['datetime'].to_list()
    cl_interval_starts = cl_interval_starts + interval_list
    
##load in relevant OMNI data. Just one year this time, but for other years will require two.

om_csv_1 = "/Users/apx059/Documents/OMNI_Raw_Data_24Yr/CSVs/2021OMNI_Raw.csv"
om_csv_2 = "/Users/apx059/Documents/OMNI_Raw_Data_24Yr/CSVs/2022OMNI_Raw.csv"

omni_1 = pd.read_csv(om_csv_1)
omni_2 = pd.read_csv(om_csv_2)

om_dfs = [omni_1, omni_2]
    
omni_all = pd.concat(om_dfs)

omni_all['datetime'] = pd.to_datetime(omni_all['datetime'],format='mixed')

omni_all = omni_all.set_index('datetime')
    
window_list = cl_interval_starts

omni_df = omni_seg(omni_all, window_list)
omni_df = omni_df.set_index('datetime')

#for df in omni_ave list:
CSV_path = '/Users/apx059/Documents/OM_SEG_TEST/'

if omni_df.size > 0:
    firstwin = omni_df.index[0]
    firstwin = str(firstwin)
    fpath = CSV_path + firstwin + 'OMNI.csv'
    omni_df.to_csv(fpath)