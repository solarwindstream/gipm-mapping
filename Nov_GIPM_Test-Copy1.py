#Apocrita script for GIPM conversion. Will require virtual env for cdflib. Adjusted for new lists
import datetime as dt
import cdflib
import pandas as pd
import numpy as np
import math

from window_det import window_det
from omni_seg_centered import omni_seg
from gipm_transform_coeffs import gipm_transform_coeffs
from GIPM_loc_conv import gipm_loc_transform
from new_xyz import new_xyz
from cone_angle_dfs import cone_angle_dfs
from math import pi

#four intervals that we're interested in first!

csv_list = ['/Users/apx059/Documents/1_Yr_Data/52_Weeks_CSVs/March-April/2001-03-29 03:33:00.027000C3.csv','/Users/apx059/Documents/1_Yr_Data/52_Weeks_CSVs/March-April/2001-03-19 14:59:00.032000C1.csv','/Users/apx059/Documents/1_Yr_Data/52_Weeks_CSVs/Dec-Feb/2002-01-17 06:41:21.192000C1.csv','/Users/apx059/Documents/1_Yr_Data/52_Weeks_CSVs/March-April/2001-03-31 12:38:00.009000C3.csv']
om_yr = pd.read_csv('/Users/apx059/Documents/1_Yr_Data/OMNI_Raw_CSVs/omni_hros_1min_20010201000000_20020201000000.csv')
om_yr['datetime'] = pd.to_datetime(om_yr['datetime'],format='mixed')
om_yr = om_yr.set_index('datetime')

def csvconv_gipm(list_of_csvs, om_yr):
    
    df_list_c1 = []
    df_list_c3 = []
    
    for i in list_of_csvs:
        if 'C1.csv' in i:
            cl_df = pd.read_csv(i)
            cl_df['datetime']=pd.to_datetime(cl_df['datetime'],format='mixed')
            cl_df = cl_df.set_index('datetime')
            df_list_c1.append(cl_df)
    
    for i in list_of_csvs:
        if 'C3.csv' in i:
            cl_df = pd.read_csv(i)
            cl_df['datetime']=pd.to_datetime(cl_df['datetime'],format='mixed')
            cl_df = cl_df.set_index('datetime')
            df_list_c3.append(cl_df)

    #determine full windows lists
    f_winds_all_1 = []

    for df in df_list_c1:
        f_winds = window_det(df)
        f_winds_all_1.append(f_winds)
    
    om_ave_list_1 = []
    
    for fw in f_winds_all_1:
        om_averages = omni_seg(om_yr, fw)
        om_averages['datetime']=pd.to_datetime(om_averages['datetime'])
        om_averages = om_averages.set_index('datetime')
        om_ave_list_1.append(om_averages)

    #find GIPM rotation matrices and scaling coefficient for every Cluster location
    GIPM_X_vec_list_1 = []
    GIPM_Y_vec_list_1 = []
    GIPM_Z_vec_list_1 = []
    FAC_coeff_list_1 = []

    for om_df in om_ave_list_1:
        GIPM_X_Vecs, GIPM_Y_Vecs, GIPM_Z_Vecs, FAC_coeffs = gipm_transform_coeffs(om_df)
        GIPM_X_vec_list_1.append(GIPM_X_Vecs)
        GIPM_Y_vec_list_1.append(GIPM_Y_Vecs)
        GIPM_Z_vec_list_1.append(GIPM_Z_Vecs)
        FAC_coeff_list_1.append(FAC_coeffs)
        
    ##print shapes of lists
    
    #for i,j,k in zip(GIPM_X_vec_list_1[0], GIPM_Y_vec_list_1[0], GIPM_Z_vec_list_1[0]):
        #print (i.shape, j.shape, k.shape)
    #for i,j,k in zip(GIPM_X_vec_list_1[1], GIPM_Y_vec_list_1[1], GIPM_Z_vec_list_1[1]):
        #print (i.shape, j.shape, k.shape)

    Cluster_GIPM_locs_list_1 = []

    for p,q,i,j,k,m in zip(f_winds_all_1, df_list_c1, GIPM_X_vec_list_1, GIPM_Y_vec_list_1, GIPM_Z_vec_list_1, FAC_coeff_list_1):
        Cluster_dt_loc = gipm_loc_transform(p,q,i,j,k,m)
        Cluster_GIPM_locs_list_1.append(Cluster_dt_loc)
        
    gipm_df_c1 = []

    for df_c1 in Cluster_GIPM_locs_list_1:
        df_c1['datetime'] = pd.to_datetime(df_c1['datetime'])
        df_c1 = df_c1.set_index('datetime')
        gipm_df_c1.append(df_c1)

    #this bit is for the new dfs 

    list_expanded_dfs_1 = []

    time_window = dt.timedelta(seconds=120)

    for i,j,k in zip(df_list_c1, gipm_df_c1, om_ave_list_1):

        cl_min_list = []
        cl_max_list = []
        cl_mean_list = []
        cone_angle_list = []
        ma_list = []
        times = []

        for m in j.index:
            start_time = m
            end_time = m + time_window
            mask = i.loc[(i.index >= start_time) & (i.index < end_time)]
            Cluster_list = mask['B_mag'].tolist()
            omni_ave_B = k.loc[m, 'B_mag']
            omni_ave_c_a = k.loc[m, 'cone angle']
            omni_ave_m_a = k.loc[m, 'M_A']
            
            if Cluster_list:
                Cluster_min = min(Cluster_list)
                Cluster_max = max(Cluster_list)
                Cluster_mean = sum(Cluster_list)/len(Cluster_list)
                Cluster_min_rat = Cluster_min/omni_ave_B
                Cluster_mean_rat = Cluster_mean/omni_ave_B
                Cluster_max_rat = Cluster_max/omni_ave_B
                
                cl_min_list.append(Cluster_min_rat)
                cl_mean_list.append(Cluster_mean_rat)
                cl_max_list.append(Cluster_max_rat)
                cone_angle_list.append(omni_ave_c_a)
                ma_list.append(omni_ave_m_a)
                times.append(m)

        if cl_min_list:
            B_val_df = pd.DataFrame({'datetime': times,'B min': cl_min_list, 'B mean': cl_mean_list, 'B max': cl_max_list, 'cone angle': cone_angle_list, 'M_A': ma_list})
            B_val_df = B_val_df.set_index('datetime')
            new_cl_df = j.join([B_val_df])
            list_expanded_dfs_1.append(new_cl_df)

        
    #for df in Cluster_GIPM_locs_list
    CSV_path = '/Users/apx059/Documents/1_Yr_Data/Nov 24 Corrected GIPM Results/'

    for df in list_expanded_dfs_1:
        if df.size > 0:
            firstwin = df.index[0]
            firstwin = str(firstwin)
            fpath = CSV_path + firstwin + 'C1.csv'
            df.to_csv(fpath)
            
    #for df in omni_ave list:

    for om_df in om_ave_list_1:
        if om_df.size > 0:
            firstwin = om_df.index[0]
            firstwin = str(firstwin)
            fpath = CSV_path + firstwin + 'OMNI_C1.csv'
            om_df.to_csv(fpath)
    
    #####C3
    
    #determine full windows lists
    f_winds_all_3 = []

    for df in df_list_c3:
        f_winds = window_det(df)
        f_winds_all_3.append(f_winds)
    
    om_ave_list_3 = []
    
    for fw in f_winds_all_3:
        om_averages = omni_seg(om_yr, fw)
        om_averages['datetime']=pd.to_datetime(om_averages['datetime'])
        om_averages = om_averages.set_index('datetime')
        om_ave_list_3.append(om_averages)

    #find GIPM rotation matrices and scaling coefficient for every Cluster location
    GIPM_X_vec_list_3 = []
    GIPM_Y_vec_list_3 = []
    GIPM_Z_vec_list_3 = []
    FAC_coeff_list_3 = []

    for om_df in om_ave_list_3:
        GIPM_X_Vecs, GIPM_Y_Vecs, GIPM_Z_Vecs, FAC_coeffs = gipm_transform_coeffs(om_df)
        GIPM_X_vec_list_3.append(GIPM_X_Vecs)
        GIPM_Y_vec_list_3.append(GIPM_Y_Vecs)
        GIPM_Z_vec_list_3.append(GIPM_Z_Vecs)
        FAC_coeff_list_3.append(FAC_coeffs)

    Cluster_GIPM_locs_list_3 = []

    for p,q,i,j,k,m in zip(f_winds_all_3, df_list_c3, GIPM_X_vec_list_3, GIPM_Y_vec_list_3, GIPM_Z_vec_list_3, FAC_coeff_list_3):
        Cluster_dt_loc = gipm_loc_transform(p,q,i,j,k,m)
        Cluster_GIPM_locs_list_3.append(Cluster_dt_loc)
        
    gipm_df_c3 = []

    for df_c3 in Cluster_GIPM_locs_list_3:
        df_c3['datetime'] = pd.to_datetime(df_c3['datetime'])
        df_c3 = df_c3.set_index('datetime')
        gipm_df_c3.append(df_c3)

    #this bit is for the new dfs 

    list_expanded_dfs_3 = []

    for i,j,k in zip(df_list_c3, gipm_df_c3, om_ave_list_3):

        cl_min_list = []
        cl_max_list = []
        cl_mean_list = []
        cone_angle_list = []
        ma_list = []
        times = []

        for m in j.index:
            start_time = m
            end_time = m + time_window
            mask = i.loc[(i.index >= start_time) & (i.index < end_time)]
            Cluster_list = mask['B_mag'].tolist()
            omni_ave_B = k.loc[m, 'B_mag']
            omni_ave_c_a = k.loc[m, 'cone angle']
            omni_ave_m_a = k.loc[m, 'M_A']
            
            if Cluster_list:
                Cluster_min = min(Cluster_list)
                Cluster_max = max(Cluster_list)
                Cluster_mean = sum(Cluster_list)/len(Cluster_list)
                Cluster_min_rat = Cluster_min/omni_ave_B
                Cluster_mean_rat = Cluster_mean/omni_ave_B
                Cluster_max_rat = Cluster_max/omni_ave_B
                
                cone_angle_list.append(omni_ave_c_a)
                ma_list.append(omni_ave_m_a)
                cl_min_list.append(Cluster_min_rat)
                cl_mean_list.append(Cluster_mean_rat)
                cl_max_list.append(Cluster_max_rat)
                
                times.append(m)

        if cl_min_list:
            B_val_df = pd.DataFrame({'datetime': times,'B min': cl_min_list, 'B mean': cl_mean_list, 'B max': cl_max_list, 'cone angle': cone_angle_list, 'M_A': ma_list})
            B_val_df = B_val_df.set_index('datetime')
            new_cl_df = j.join([B_val_df])
            list_expanded_dfs_3.append(new_cl_df)

    for df in list_expanded_dfs_3:
        if df.size > 0:
            firstwin = df.index[0]
            firstwin = str(firstwin)
            fpath = CSV_path + firstwin + 'C3.csv'
            df.to_csv(fpath)

    for om_df in om_ave_list_3:
        if om_df.size > 0:
            firstwin = om_df.index[0]
            firstwin = str(firstwin)
            fpath = CSV_path + firstwin + 'OMNI_C3.csv'
            om_df.to_csv(fpath)
        
    #########################
    
csvconv_gipm(csv_list, om_yr)

