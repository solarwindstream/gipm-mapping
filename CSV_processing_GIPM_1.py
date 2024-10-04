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
from cone_angle_dfs import cone_angle_dfs

#run a test of the GIPM conversion for six weeks of CSVs

#open cluster CSVs
path = r'/data/scratch/apx059/Feb-March-2001/**'
#path = r'/data/scratch/apx059/48_Wks_CSVs/Dec-Feb/**'
#path = r'/data/scratch/apx059/48_Wks_CSVs/March-April/**'
#path = r'/data/scratch/apx059/48_Wks_CSVs/May-June/**'
#path = r'/data/scratch/apx059/48_Wks_CSVs/July-November/**'

list_all = []
for path1 in glob.glob(path, recursive=True):
    list_all.append(path1)

#list with only files, not folders

cluster_csv_list = []

for element in list_all:
    if '.csv' in element:
        cluster_csv_list.append(element)

cluster_dfs = []

for file in cluster_csv_list:
    df = pd.read_csv(file)
    df['datetime'] = pd.to_datetime(df['datetime'])
    df = df.set_index('datetime')
    cluster_dfs.append(df)

#append OMNI dfs to each other
path = r'/data/scratch/apx059/OMNI_Raw/**'

list_all = []
for path1 in glob.glob(path, recursive=True):
    list_all.append(path1)

#list with only files, not folders

omni_csv_list = []

for element in list_all:
    if '.csv' in element:
        omni_csv_list.append(element)

omni_df_list = []

for file in omni_csv_list:
    df = pd.read_csv(file)
    df['datetime'] = pd.to_datetime(df['datetime'])
    df = df.set_index('datetime')
    omni_df_list.append(df)

omni_df = pd.concat(omni_df_list)

#omni_df = pd.read_csv(omni_csv)
#omni_df['datetime'] = pd.to_datetime(omni_df['datetime'])
#omni_df = omni_df.set_index('datetime')


#determine full windows lists
f_winds_all = []

for df in cluster_dfs:
    f_winds = window_det(df)
    f_winds_all.append(f_winds)

om_ave_dfs_list = []

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
    
#split by cone angle
rad_df_list = []
spir_df_list = []
perp_df_list = []

for om_ave, cl_df, GIPM_lc in zip(om_ave_dfs_list, cluster_dfs, Cluster_GIPM_loc_all):
    rad_df, spir_df, perp_df = cone_angle_dfs(om_ave,cl_df, GIPM_lc)
    rad_df_list.append(rad_df)
    spir_df_list.append(spir_df)
    perp_df_list.append(perp_df)
    
#save files as CSVs

#omni first. add numbers to allow saved dupes for multi sc

fw_list = []

for df in om_ave_dfs_list:
    firstwin = df.index[0]
    firstwin = str(firstwin)
    if firstwin not in fw_list:
        fpath = '/data/scratch/apx059/48_Wks_Res/omni_ave_' + firstwin + '.csv'
        df.to_csv(fpath)
    else:
        counter = fw_list.count(firstwin)
        counter = str(counter)
        fpath = '/data/scratch/apx059/48_Wks_Res/omni_ave_' + firstwin + counter + '.csv'
        df.to_csv(fpath) 
    fw_list.append(firstwin)
    
#now Cluster
fw_list = []

for df in rad_df_list:
    if df.size > 0:
        firstwin = df.loc[0, 'window start']
        firstwin = str(firstwin)
        if firstwin not in fw_list:
            fpath = '/data/scratch/apx059/48_Wks_Res/cluster_rad' + firstwin + '.csv'
            df.to_csv(fpath)
        else:
            counter = fw_list.count(firstwin)
            counter = str(counter)
            fpath = '/data/scratch/apx059/48_Wks_Res/cluster_rad' + firstwin + counter + '.csv'
            df.to_csv(fpath) 
    fw_list.append(firstwin)

fw_list = []
    
for df in spir_df_list:
    if df.size > 0:
        firstwin = df.loc[0, 'window start']
        firstwin = str(firstwin)
        if firstwin not in fw_list:
            fpath = '/data/scratch/apx059/48_Wks_Res/cluster_spir' + firstwin + '.csv'
            df.to_csv(fpath)
        else:
            counter = fw_list.count(firstwin)
            counter = str(counter)
            fpath = '/data/scratch/apx059/48_Wks_Res/cluster_spir' + firstwin + counter + '.csv'
            df.to_csv(fpath) 
    fw_list.append(firstwin)

fw_list = []

for df in perp_df_list:
    if df.size > 0:
        firstwin = df.loc[0, 'window start']
        firstwin = str(firstwin)
        if firstwin not in fw_list:
            fpath = '/data/scratch/apx059/48_Wks_Res/cluster_perp' + firstwin + '.csv'
            df.to_csv(fpath)
        else:
            counter = fw_list.count(firstwin)
            counter = str(counter)
            fpath = '/data/scratch/apx059/48_Wks_Res/cluster_perp' + firstwin + counter + '.csv'
            df.to_csv(fpath) 
    fw_list.append(firstwin)
    