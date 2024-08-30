#for use on cluster
import pandas as pd
import numpy as np
import datetime as dt
from omni_seg import omni_seg
from gipm_transform_coeffs import gipm_transform_coeffs
from new_xyz import new_xyz
import glob
from math import pi
import math
from gipm_locs_quick import gipm_locs_quick

#open OMNI CSVs
list_all = []

path = r"/data/scratch/apx059/OMNI_Raw/**"

for path in glob.glob(path, recursive=True):
    list_all.append(path)
    
#list with only files, not folders
om_csvs = []

for element in list_all:
    if '.csv' in element:
        om_csvs.append(element)
        
om_dfs = []

for element in om_csvs:
    om = pd.read_csv(element)
    om_dfs.append(om)
    
omni_all = pd.concat(om_dfs)
omni_all['datetime'] = pd.to_datetime(omni_all['datetime'])

omni_all = omni_all.set_index('datetime')

for element in om_dfs:
    element['datetime'] = pd.to_datetime(element['datetime'])
    element = element.set_index('datetime', inplace = True)

##load Cluster CSVs

list_all = []

path = r'/data/scratch/apx059/Cluster_Raw/**'

for path in glob.glob(path, recursive=True):
    list_all.append(path)
    
#list with only files, not folders
cl_file_list = []

for element in list_all:
    if '.csv' in element:
        cl_file_list.append(element)

cl_dfs = []

for file in cl_file_list:
    df = pd.read_csv(file,encoding='utf-8')
    df['datetime'] = pd.to_datetime(df['datetime'])
    df.set_index('datetime', inplace = True)
    cl_dfs.append(df)
    
#need to reorder bc they're not automatically in the right order
reordered_om_df = []

for df in cl_dfs:
    fd_list = df.index
    fd = fd_list[0]
    first_read = fd.strftime("%Y")
    for omni_temp in om_dfs:
        fd_om_list_temp = omni_temp.index
        fd_om_temp = fd_om_list_temp[0]
        first_read_om_t = fd_om_temp.strftime("%Y")
        if first_read_om_t == first_read:
            reordered_om_df.append(omni_temp)
            
#break into intervals with averaged OMNI values. for subsequent GIPM plotting
#first, pair Cluster and OMNI dataframes

om_ave_list = []
dayside_cl_list = []

for df, om_df in zip(cl_dfs,reordered_om_df):
    #mask to just dayside
    cluster_df_day = df.loc[(df['X_gse']> 0)]
    datetime_list = cluster_df_day.index
    om_averages = omni_seg(om_df, datetime_list)
    om_ave_list.append(om_averages)
    dayside_cl_list.append(cluster_df_day)
    
for om_ave_df in om_ave_list:
    CSV_path = '/data/scratch/apx059/Omni_Aves_SSC'
    start_ind = om_ave_df.index[0]
    start_year = start_ind.strftime("%Y")
    start_year_int = int(start_year)
    end_year_int = start_year_int + 1
    end_year = str(end_year_int)
    CSV_file_name = CSV_path + 'OMNI' + '_Feb' + start_year + '_1yr' + '.csv'
    om_ave_df.to_csv(CSV_file_name)

df_no = 0
for cl_day in dayside_cl_list:
    CSV_path = '/data/scratch/apx059/Cluster_Dayside_SSC'
    df_ref = str(df_no)
    CSV_file_name = CSV_path + 'Cluster Dayside' + '_dfref' + df_ref +'.csv'
    cl_day.to_csv(CSV_file_name)
    df_no = df_no + 1

#find GIPM rotation matrices and scaling coefficient for every Cluster location
GIPM_mat_list = []
FAC_coeff_list = []

for om_df in om_ave_list:
    GIPM_mat, FAC_coeffs = gipm_transform_coeffs(om_df)
    GIPM_mat_list.append(GIPM_mat)
    FAC_coeff_list.append(FAC_coeffs)
    
Cluster_GIPM_locs_list = []

for i,j,k in zip(dayside_cl_list, GIPM_mat_list, FAC_coeff_list):
    Cluster_dt_loc = gipm_locs_quick(i, j, k)
    Cluster_GIPM_locs_list.append(Cluster_dt_loc)
    
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