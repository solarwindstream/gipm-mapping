import csv
import pandas as pd
import numpy as np
import glob
import datetime as dt
import math

from window_det import window_det
from omni_seg_centered import omni_seg
from gipm_transform_coeffs import gipm_transform_coeffs
from GIPM_loc_conv import gipm_loc_transform

#run a test of the GIPM conversion for six weeks of CSVs

#open cluster CSVs from Scratch storage
path = r'/data/scratch/apx059/1_Yr_Data/Cluster_CSVs_2001-02-01-2001-03-15'

list_all = []
for path1 in glob.glob(path, recursive=True):
    list_all.append(path1)

#list with only files, not folders
cluster_csv_list = []

for element in list_all:
    if '.csv' in element:
        cluster_csv_list.append(element)

for file in cluster_csv_list:
    df = pd.read_csv(file)
    df['datetime'] = pd.to_datetime(df['datetime'])
    df = df.set_index('datetime')
    cluster_dfs.append(df)

#append OMNI dfs to each other

omni_csv = r'/data/scratch/apx059/OMNI_Raw_CSVs/omni_hros_1min_20010201000000_20020201000000.csv'
#omni_df_list = []

#for file in omni_csv_list:
#    df = pd.read_csv(file)
#    df['datetime'] = pd.to_datetime(df['datetime'])
#    df = df.set_index('datetime')
#    omni_df_list.append(df)

#omni_df = pd.concat(omni_df_list)

omni_df = pd.read_csv(omni_csv)
omni_df['datetime'] = pd.to_datetime(omni_df['datetime'])
omni_df = omni_df.set_index('datetime')


#determine full windows lists
f_winds_all = []

for df in cluster_dfs:
    f_winds = window_det(df)
    f_winds_all.append(f_winds)

for wind_list in f_winds_all:
    om_ave = omni_seg(omni_df, wind_list)
    om_ave['datetime']=pd.to_datetime(om_ave['datetime'])
    om_ave = om_ave.set_index('datetime')
    om_ave_dfs_list.append(om_ave)
    
#find GIPM rotation matrices and scaling coefficient for every two minute interval
GIPM_mat_list = []
FAC_list = []

for df in om_ave_dfs_list:
    GIPM_mat, FAC_coeffs = gipm_transform_coeffs(df)
    GIPM_mat_list.append(GIPM_mat)
    FAC_list.append(FAC_coeffs)
    
#use GIPM matrices and scaling coeff to generate average GIPM location for spacecraft in 
#each two minute interval.

Cluster_GIPM_loc_all = []

for f_w,cl_df,GIPM_mat,FAC_c in zip(f_winds_all, cluster_dfs, GIPM_mat_list, FAC_list):
    Cluster_GIPM_loc = gipm_loc_transform(f_w, cl_df, GIPM_mat, FAC_c)
    Cluster_GIPM_loc_all.append(Cluster_GIPM_loc)
    
    