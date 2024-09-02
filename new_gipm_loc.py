#updated version of script for 23 yr processing

import pandas as pd
import numpy as np
import datetime as dt
from omni_seg import omni_seg
from gipm_transform_coeffs import gipm_transform_coeffs
from new_xyz import new_xyz
import glob
from math import pi
import math
from XMA_finder import XMA_finder
from histo_plot import histo_plot
from gipm_locs_quick import gipm_locs_quick

#load dayside cl list and om ave list

list_all = []

path = "/data/scratch/apx059/Cluster_Dayside_SSC/**"

for path in glob.glob(path, recursive=True):
    list_all.append(path)
    
#list with only files, not folders
cl_ds_csvs = []

for element in list_all:
    if '.csv' in element:
        cl_ds_csvs.append(element)
        
cl_ds_dfs = []

for element in cl_ds_csvs:
    clu = pd.read_csv(element)
    cl_ds_dfs.append(clu)

for element in cl_ds_dfs:
    element['datetime'] = pd.to_datetime(element['datetime'])
    element = element.set_index('datetime', inplace = True)
    
######################

list_all = []

path = "/data/scratch/apx059/Omni_Aves_SSC/**"

for path in glob.glob(path, recursive=True):
    list_all.append(path)
    
#list with only files, not folders
om_csvs = []

for element in list_all:
    if '.csv' in element:
        om_csvs.append(element)
        
om_ave_dfs = []

for element in om_csvs:
    om = pd.read_csv(element)
    om_ave_dfs.append(om)

for element in om_ave_dfs:
    element['datetime'] = pd.to_datetime(element['datetime'])
    element = element.set_index('datetime', inplace = True)

#first identify WHICH of the omni average dfs it matches up with, find the list index for it, 
#then find the corresponding mat & FAC entries and do transform.

#list by first year in OMNI dfs
omni_years = []

for omni_df in om_ave_dfs:
    fd_om_list_temp = omni_df.index
    fd_om_temp = fd_om_list_temp[0]
    first_read_om_t = fd_om_temp.strftime("%Y")
    omni_years.append(first_read_om_t)

print('omni years')

ind_no_list = []

for df in cl_ds_dfs:
    fd_list = df.index
    fd = fd_list[0]
    first_read = fd.strftime("%Y")
    ind_no = omni_years.index(first_read)
    ind_no_list.append(ind_no)
    
Cluster_GIPM_locs_list = []

for i,j in zip(cl_ds_dfs, ind_no_list):
    mat_inst = GIPM_mat_list[j]
    FAC_inst = FAC_coeff_list[j]
    Cluster_dt_loc = gipm_locs_quick(i, mat_inst, FAC_inst)
    Cluster_GIPM_locs_list.append(Cluster_dt_loc)
    print('done!')
    
#for df in Cluster_GIPM_locs_list
CSV_path = '/data/scratch/apx059/Cluster_GIPM_SSC/'
startref = 1
array_list = []

for df in Cluster_GIPM_locs_list:
    #formatting numpy array correctly
    list_loc = df['GIPM Loc'].tolist()
    array_loc = np.asarray(list_loc)
    list_time = df['datetime'].tolist()
    array_time = np.asarray(list_time)
    array_time = [[item] for item in array_time]
    array_both = np.append(array_time, array_loc, axis=1)
    array_list.append(array_both)
    #saving file
    ref = str(startref)
    start_year = list_time[0].strftime("%Y")
    start_year_int = int(start_year)
    end_year_int = start_year_int + 1
    end_year = str(end_year_int)
    filename = CSV_path + '_Feb' + start_year + '_Feb' + end_year + 'df_' + ref + '.csv'
    np.savetxt(filename,array_both, fmt='%s,%f,%f,%f')
    startref = startref + 1