#Apocrita script for GIPM conversion. Will require virtual env for cdflib. Adjusted for new lists
import datetime as dt
import cdflib
import pandas as pd
import numpy as np
import math

from C1_Cluster_CDF_conv import C1_cdf_conv
from C2_Cluster_CDF_conv import C2_cdf_conv
from C3_Cluster_CDF_conv import C3_cdf_conv
from C4_Cluster_CDF_conv import C4_cdf_conv

from window_det import window_det
from omni_seg_centered import omni_seg
from gipm_transform_coeffs import gipm_transform_coeffs
from GIPM_loc_conv import gipm_loc_transform
from new_xyz import new_xyz
from cone_angle_dfs import cone_angle_dfs
from math import pi

cdf_list_path = '/data/scratch/apx059/23_Years_Data/CDFs/CDFs_'

def cdfconv_gipm(path, year):
    
    #now make those cdfs into CSVs!
    list_files = pd.read_csv(path, header=None)
    #list_files_T = list_files.T
    #print('list_files_T:', list_files_T)
    #list_files_mask = list_files.loc[((list_files.index >= batch_start) & (list_files.index < batch_end))]
    list_files = list_files.rename(columns={0:'fnames'})
    list_cdfs = list_files['fnames'].to_list()
    year_path = '/data/scratch/apx059/23_Years_Data/CDFs/CDFs_' + year_n + '/'
    
    #now iterate over list to convert

    df_list_c1 = []
    df_list_c2 = []
    df_list_c3 = []
    df_list_c4 = []

    for i in list_cdfs:
        if 'C1' in i:
            fpath = year_path + i
            df_c1 = C1_cdf_conv(fpath)
            a = df_c1.empty
            if not a:
                df_list_c1.append(df_c1)
        if 'C2' in i:
            fpath = year_path + i
            df_c2 = C2_cdf_conv(fpath)
            a = df_c2.empty
            if not a:
                df_list_c2.append(df_c2)
        if 'C3' in i:   
            fpath = year_path + i
            df_c3 = C3_cdf_conv(fpath)
            a = df_c3.empty
            if not a:
                df_list_c3.append(df_c3) 
        if 'C4' in i:   
            fpath = year_path + i
            df_c4 = C4_cdf_conv(fpath)
            a = df_c4.empty
            if not a:
                df_list_c4.append(df_c4)

    #list all OMNI files, but only import relevant ones
    #load ONLY this year's omni df
    
    omni_csv_list_path = r'/data/home/apx059/gipm-mapping/listofomnis.csv'
    omni_files = pd.read_csv(omni_csv_list_path, header=None)
    omni_files = omni_files.rename(columns={0:'fnames'})
    list_omni_csvs = omni_files['fnames'].to_list()
    #res doesn't work because omni date time strings contain non-date repeats
    #res = [i for i in list_omni_csvs if year_n in i]
    #omni_hros_1min_20200201000000_20210201000000.csv
    #omni_hros_1min_20190201000000_20200201000000.csv
    omni_path_1 = '/data/scratch/apx059/OMNI_Raw/' + 'omni_hros_1min_20020201000000_20030201000000.csv'
    omni_path_2 = '/data/scratch/apx059/OMNI_Raw/' + 'omni_hros_1min_20030201000000_20040201000000.csv'
    om_1 = pd.read_csv(omni_path_1)
    om_2 = pd.read_csv(omni_path_2)
    om = pd.concat([om_1, om_2])

    om['datetime'] = pd.to_datetime(om['datetime'])
    om = om.set_index('datetime')

    ##########################
    #GIPM conversion!
    ####PRESERVE SC LISTINGS
    ####SC1

    #determine full windows lists
    f_winds_all_1 = []

    for df in df_list_c1:
        f_winds = window_det(df)
        f_winds_all_1.append(f_winds)
    
    om_ave_list_1 = []
    
    for fw in f_winds_all_1:
        om_averages = omni_seg(om, fw)
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
        GIPM_mat_list_1.append(GIPM_mat)
        FAC_coeff_list_1.append(FAC_coeffs)

    Cluster_GIPM_locs_list_1 = []

    for p,i,j,k in zip(f_winds_all_1, df_list_c1, GIPM_mat_list_1, FAC_coeff_list_1):
        Cluster_dt_loc = gipm_loc_transform(p, i, j, k)
        Cluster_GIPM_locs_list_1.append(Cluster_dt_loc)
        
    #for df in Cluster_GIPM_locs_list
    CSV_path = '/data/scratch/apx059/23_Years_Data/CSVs/Updated_GIPM/'

    for df in Cluster_GIPM_locs_list_1:
        if df.size > 0:
            firstwin = df.loc[0, 'datetime']
            firstwin = str(firstwin)
            fpath = CSV_path + firstwin + 'C1.csv'
            df.to_csv(fpath)
                
    #for df in omni_ave list:
    CSV_path = '/data/scratch/apx059/23_Years_Data/CSVs/OMNI_Aves/'

    for om_df in om_ave_list_1:
        if om_df.size > 0:
            firstwin = om_df.index[0]
            firstwin = str(firstwin)
            fpath = CSV_path + firstwin + 'OMNI_C1.csv'
            om_df.to_csv(fpath)

    #####C2!!!!
    #determine full windows lists
    f_winds_all_2 = []

    for df in df_list_c2:
        f_winds = window_det(df)
        f_winds_all_2.append(f_winds)
    
    om_ave_list_2 = []
    
    for fw in f_winds_all_2:
        om_averages = omni_seg(om, fw)
        om_averages['datetime']=pd.to_datetime(om_averages['datetime'])
        om_averages = om_averages.set_index('datetime')
        om_ave_list_2.append(om_averages)

    #find GIPM rotation matrices and scaling coefficient for every Cluster location
    GIPM_mat_list_2 = []
    FAC_coeff_list_2 = []

    for om_df in om_ave_list_2:
        GIPM_mat, FAC_coeffs = gipm_transform_coeffs(om_df)
        GIPM_mat_list_2.append(GIPM_mat)
        FAC_coeff_list_2.append(FAC_coeffs)

    Cluster_GIPM_locs_list_2 = []

    for p,i,j,k in zip(f_winds_all_2, df_list_c2, GIPM_mat_list_2, FAC_coeff_list_2):
        Cluster_dt_loc = gipm_loc_transform(p, i, j, k)
        Cluster_GIPM_locs_list_2.append(Cluster_dt_loc)
        
    #for df in Cluster_GIPM_locs_list
    CSV_path = '/data/scratch/apx059/23_Years_Data/CSVs/Updated_GIPM/'
    array_list = []

    for df in Cluster_GIPM_locs_list_2:
        if df.size > 0:
            firstwin = df.loc[0, 'datetime']
            firstwin = str(firstwin)
            fpath = CSV_path + firstwin + 'C2.csv'
            df.to_csv(fpath)
    
    #for df in omni_ave list:
    CSV_path = '/data/scratch/apx059/23_Years_Data/CSVs/OMNI_Aves/'
           
    for om_df in om_ave_list_2:
        if om_df.size > 0:
            firstwin = om_df.index[0]
            firstwin = str(firstwin)
            fpath = CSV_path + firstwin + 'OMNI_C2.csv'
            om_df.to_csv(fpath)
                
    #####C3
    #determine full windows lists
    f_winds_all_3 = []

    for df in df_list_c3:
        f_winds = window_det(df)
        f_winds_all_3.append(f_winds)
    
    om_ave_list_3 = []
    
    for fw in f_winds_all_3:
        om_averages = omni_seg(om, fw)
        om_averages['datetime']=pd.to_datetime(om_averages['datetime'])
        om_averages = om_averages.set_index('datetime')
        om_ave_list_3.append(om_averages)

    #find GIPM rotation matrices and scaling coefficient for every Cluster location
    GIPM_mat_list_3 = []
    FAC_coeff_list_3 = []

    for om_df in om_ave_list_3:
        GIPM_mat, FAC_coeffs = gipm_transform_coeffs(om_df)
        GIPM_mat_list_3.append(GIPM_mat)
        FAC_coeff_list_3.append(FAC_coeffs)

    Cluster_GIPM_locs_list_3 = []

    for p,i,j,k in zip(f_winds_all_3, df_list_c3, GIPM_mat_list_3, FAC_coeff_list_3):
        Cluster_dt_loc = gipm_loc_transform(p, i, j, k)
        Cluster_GIPM_locs_list_3.append(Cluster_dt_loc)
        
    #for df in Cluster_GIPM_locs_list
    CSV_path = '/data/scratch/apx059/23_Years_Data/CSVs/Updated_GIPM/'
    array_list = []

    for df in Cluster_GIPM_locs_list_3:
        if df.size > 0:
            firstwin = df.loc[0, 'datetime']
            firstwin = str(firstwin)
            fpath = CSV_path + firstwin + 'C3.csv'
            df.to_csv(fpath)
    
    #for df in omni_ave list:
    CSV_path = '/data/scratch/apx059/23_Years_Data/CSVs/OMNI_Aves/'

    for om_df in om_ave_list_3:
        if om_df.size > 0:
            firstwin = om_df.index[0]
            firstwin = str(firstwin)
            fpath = CSV_path + firstwin + 'OMNI_C3.csv'
            om_df.to_csv(fpath)
        
    #####C4!!!!
    #determine full windows lists
    f_winds_all_4 = []

    for df in df_list_c4:
        f_winds = window_det(df)
        f_winds_all_4.append(f_winds)
    
    om_ave_list_4 = []
    
    for fw in f_winds_all_4:
        om_averages = omni_seg(om, fw)
        om_averages['datetime']=pd.to_datetime(om_averages['datetime'])
        om_averages = om_averages.set_index('datetime')
        om_ave_list_4.append(om_averages)

    #find GIPM rotation matrices and scaling coefficient for every Cluster location
    GIPM_mat_list_4 = []
    FAC_coeff_list_4 = []

    for om_df in om_ave_list_4:
        GIPM_mat, FAC_coeffs = gipm_transform_coeffs(om_df)
        GIPM_mat_list_4.append(GIPM_mat)
        FAC_coeff_list_4.append(FAC_coeffs)

    Cluster_GIPM_locs_list_4 = []

    for p,i,j,k in zip(f_winds_all_4, df_list_c4, GIPM_mat_list_4, FAC_coeff_list_4):
        Cluster_dt_loc = gipm_loc_transform(p, i, j, k)
        Cluster_GIPM_locs_list_4.append(Cluster_dt_loc)
        
    #for df in Cluster_GIPM_locs_list
    CSV_path = '/data/scratch/apx059/23_Years_Data/CSVs/Updated_GIPM/'
    array_list = []

    for df in Cluster_GIPM_locs_list_4:
        if df.size > 0:
            firstwin = df.loc[0, 'datetime']
            firstwin = str(firstwin)
            fpath = CSV_path + firstwin + 'C4.csv'
            df.to_csv(fpath)
    
    #for df in omni_ave list:
    CSV_path = '/data/scratch/apx059/23_Years_Data/CSVs/OMNI_Aves/'

    for om_df in om_ave_list_4:
        if om_df.size > 0:
            firstwin = om_df.index[0]
            firstwin = str(firstwin)
            fpath = CSV_path + firstwin + 'OMNI_C4.csv'
            om_df.to_csv(fpath)
            
    #########################
    
year_n = '2003'
list_of_cdfs = cdf_list_path + year_n + '.csv'
cdfconv_gipm(list_of_cdfs, year_n)

