#Apocrita script for 4min GIPM conversion w/spectra INCLUDED. Will require virtual env for cdflib.

import datetime as dt
import cdflib
import pandas as pd
import numpy as np
import math

from Cluster_CDF_conv import Cluster_cdf_conv

from window_det_4min import window_det
from omni_seg_centered import omni_seg
from gipm_transform_coeffs import gipm_transform_coeffs_mean, gipm_transform_coeffs_median
from GIPM_loc_conv import gipm_loc_transform
from new_xyz import new_xyz
from cone_angle_dfs import cone_angle_dfs
from math import pi
import statistics
from FFT_Hann import FFT_Hann
import glob

pd.options.mode.chained_assignment = None

cdf_list_path = '/data/scratch/apx059/23_Years_Data/CDFs/CDFs_'
year_n = '2003'

def cdfconv_gipm(year):
    
    #rework to use glob
    
    #folder with input CDFs

    path = '/data/SPCS-HIETALA-Shocks/GIPM-MAPPING/Cluster_CDFs/CDFs_' + year + '/Full_CDFs/**'
    
    list_cdfs = []
    
    for path in glob.glob(path, recursive=True):
        if '.cdf' in path:
            list_cdfs.append(path)
            
    #now put each filepath into conversion, checking spacecraft number.
    
    df_list_c1 = []
    df_list_c2 = []
    df_list_c3 = []
    df_list_c4 = []

    for i in list_cdfs:
        if 'C1' in i:
            df_c1 = Cluster_cdf_conv(i, 'C1')
            a = df_c1.empty
            if not a:
                df_list_c1.append(df_c1)
        if 'C2' in i:
            df_c2 = Cluster_cdf_conv(i, 'C2')
            a = df_c2.empty
            if not a:
                df_list_c2.append(df_c2)
        if 'C3' in i:   
            df_c3 = Cluster_cdf_conv(i, 'C3')
            a = df_c3.empty
            if not a:
                df_list_c3.append(df_c3) 
        if 'C4' in i:   
            df_c4 = Cluster_cdf_conv(i, 'C4')
            a = df_c4.empty
            if not a:
                df_list_c4.append(df_c4)

    #load ONLY this year's omni dfs
    
    omni_path = '/data/SPCS-HIETALA-Shocks/GIPM-MAPPING/New_OMNI_Raw_Files/CSVs/Raw_OMNI_' + year + '.csv'
    om = pd.read_csv(omni_path)

    om['datetime'] = pd.to_datetime(om['datetime'])
    om = om.set_index('datetime')

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
    
    GIPM_X_mean_list_1 = []
    GIPM_Y_mean_list_1 = []
    GIPM_Z_mean_list_1 = []
    FAC_coeff_mean_list_1 = []
    
    GIPM_X_median_list_1 = []
    GIPM_Y_median_list_1 = []
    GIPM_Z_median_list_1 = []
    FAC_coeff_median_list_1 = []

    for om_df in om_ave_list_1:
        
        #mean
        GIPM_X_Vecs, GIPM_Y_Vecs, GIPM_Z_Vecs, FAC_coeffs = gipm_transform_coeffs_mean(om_df)
        GIPM_X_mean_list_1.append(GIPM_X_Vecs)
        GIPM_Y_mean_list_1.append(GIPM_Y_Vecs)
        GIPM_Z_mean_list_1.append(GIPM_Z_Vecs)
        FAC_coeff_mean_list_1.append(FAC_coeffs)
        
        #median
        GIPM_X_Vecs, GIPM_Y_Vecs, GIPM_Z_Vecs, FAC_coeffs = gipm_transform_coeffs_median(om_df)
        GIPM_X_median_list_1.append(GIPM_X_Vecs)
        GIPM_Y_median_list_1.append(GIPM_Y_Vecs)
        GIPM_Z_median_list_1.append(GIPM_Z_Vecs)
        FAC_coeff_median_list_1.append(FAC_coeffs)

    
    #now combine these basis vectors and co-efficients with Cluster GSE location data to transform into GIPM
    Cluster_GIPM_locs_list_1 = []
    
    #mean AND median locations both produced by gipm_loc_transform
    for p,q,i,j,k,m,n,r,s,t in zip(f_winds_all_1, df_list_c1, GIPM_X_mean_list_1, GIPM_Y_mean_list_1, GIPM_Z_mean_list_1, FAC_coeff_mean_list_1, GIPM_X_median_list_1, GIPM_Y_median_list_1, GIPM_Z_median_list_1, FAC_coeff_median_list_1):
        Cluster_dt_loc = gipm_loc_transform(p,q,i,j,k,m,n,r,s,t)
        Cluster_GIPM_locs_list_1.append(Cluster_dt_loc)
        
    gipm_df_c1 = []

    for df_c1 in Cluster_GIPM_locs_list_1:
        df_c1['datetime'] = pd.to_datetime(df_c1['datetime'])
        df_c1 = df_c1.set_index('datetime')
        gipm_df_c1.append(df_c1)

    #now find min, max and mean B, integrated power, and save with relevant OMNI average values.

    list_expanded_dfs_1 = []

    time_window = dt.timedelta(seconds=240)

    for i,j,k in zip(df_list_c1, gipm_df_c1, om_ave_list_1):

        #empty lists set up for all variables
        cl_min_list, cl_max_list, cl_mean_list, cl_median_list, cl_std_list, cone_angle_list, ma_list, max_IMF_dev_list, omni_B_list, para_i_p_list, perp_i_p_list, times = ([] for i in range(12))
        
        for m in j.index:
            start_time = m
            end_time = m + time_window
            mask = i.loc[(i.index >= start_time) & (i.index < end_time)]
            Cluster_list = mask['B_mag'].tolist()
            omni_ave_B = k.loc[m, 'B_mag']
            omni_ave_c_a = k.loc[m, 'cone angle']
            omni_ave_m_a = k.loc[m, 'M_A']
            omni_ave_max_IMF_dev = k.loc[m, 'max cone angle deviation']
            
            if Cluster_list:
                Cluster_min = min(Cluster_list)
                Cluster_max = max(Cluster_list)
                Cluster_median = statistics.median(Cluster_list)
                Cluster_mean = sum(Cluster_list)/len(Cluster_list)
                Cluster_stf = statistics.stdev(Cluster_list)
                Cluster_stf_rat = Cluster_stf/omni_ave_B
                Cluster_min_rat = Cluster_min/omni_ave_B
                Cluster_mean_rat = Cluster_mean/omni_ave_B
                Cluster_max_rat = Cluster_max/omni_ave_B
                Cluster_median_rat = Cluster_median/omni_ave_B
                
                #4 minute power spectrum
                para_int_p, perp_int_p, freq, power_s_para, power_s_perp_1, power_s_perp_2 = FFT_Hann(mask)
                #save freq and power_s as their own df, named after 1- window start and 2- sc
                spectral_df = pd.DataFrame({'Freq':freq, 'Parallel Power':power_s_para, 'Perp 1 Power': power_s_perp_1, 'Perp 2 Power':power_s_perp_2})
                spec_fname = str(m)+'FS_C1.csv'
                fpath_spec = '/data/scratch/apx059/23_Years_Data/CSVs/GIPM_4mins/Fourier_Products/Fourier_Spectra_'+year_n+'/'+spec_fname
                spectral_df.to_csv(fpath_spec)
                
                cl_min_list.append(Cluster_min_rat)
                cl_mean_list.append(Cluster_mean_rat)
                cl_max_list.append(Cluster_max_rat)
                cl_median_list.append(Cluster_median_rat)
                cl_std_list.append(Cluster_stf_rat)
                cone_angle_list.append(omni_ave_c_a)
                para_i_p_list.append(para_int_p)
                perp_i_p_list.append(perp_int_p)
                ma_list.append(omni_ave_m_a)
                omni_B_list.append(omni_ave_B)
                times.append(m)

        if cl_min_list:
            B_val_df = pd.DataFrame({'datetime': times,'B min': cl_min_list, 'B mean': cl_mean_list, 'B max': cl_max_list, 'B median': cl_median_list, 'B standard deviation':cl_std_list, 'cone angle': cone_angle_list, 'M_A': ma_list, 'ULF Parallel Power':para_i_p_list, 'ULF Perpendicular Power': perp_i_p_list})
            B_val_df = B_val_df.set_index('datetime')
            new_cl_df = j.join([B_val_df])
            list_expanded_dfs_1.append(new_cl_df)

        
    #for df in Cluster_GIPM_locs_list
    CSV_path = '/data/scratch/apx059/23_Years_Data/CSVs/GIPM_4mins/'+year_n+'_Integrated_CSVs/'
    OMNI_CSV_path = '/data/scratch/apx059/23_Years_Data/CSVs/GIPM_4mins/'+year_n+'_OMNI_CSVs/'

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
            fpath = OMNI_CSV_path + firstwin + 'OMNI_C1.csv'
            om_df.to_csv(fpath)
            
    #####C2
    
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
    GIPM_X_vec_list_2 = []
    GIPM_Y_vec_list_2 = []
    GIPM_Z_vec_list_2 = []
    FAC_coeff_list_2 = []

    for om_df in om_ave_list_2:
        GIPM_X_Vecs, GIPM_Y_Vecs, GIPM_Z_Vecs, FAC_coeffs = gipm_transform_coeffs(om_df)
        GIPM_X_vec_list_2.append(GIPM_X_Vecs)
        GIPM_Y_vec_list_2.append(GIPM_Y_Vecs)
        GIPM_Z_vec_list_2.append(GIPM_Z_Vecs)
        FAC_coeff_list_2.append(FAC_coeffs)

    Cluster_GIPM_locs_list_2 = []

    for p,q,i,j,k,m in zip(f_winds_all_2, df_list_c2, GIPM_X_vec_list_2, GIPM_Y_vec_list_2, GIPM_Z_vec_list_2, FAC_coeff_list_2):
        Cluster_dt_loc = gipm_loc_transform(p,q,i,j,k,m)
        Cluster_GIPM_locs_list_2.append(Cluster_dt_loc)
        
    gipm_df_c2 = []

    for df_c2 in Cluster_GIPM_locs_list_2:
        df_c2['datetime'] = pd.to_datetime(df_c2['datetime'])
        df_c2 = df_c2.set_index('datetime')
        gipm_df_c2.append(df_c2)

    #this bit is for the new dfs 

    list_expanded_dfs_2 = []

    for i,j,k in zip(df_list_c2, gipm_df_c2, om_ave_list_2):

        cl_min_list = []
        cl_max_list = []
        cl_mean_list = []
        cl_median_list = []
        cl_std_list = []
        cone_angle_list = []
        para_i_p_list = []
        perp_i_p_list = []
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
                Cluster_median = statistics.median(Cluster_list)
                Cluster_mean = sum(Cluster_list)/len(Cluster_list)
                Cluster_stf = statistics.stdev(Cluster_list)
                Cluster_stf_rat = Cluster_stf/omni_ave_B
                Cluster_min_rat = Cluster_min/omni_ave_B
                Cluster_mean_rat = Cluster_mean/omni_ave_B
                Cluster_max_rat = Cluster_max/omni_ave_B
                Cluster_median_rat = Cluster_median/omni_ave_B
                
                #4 minute power spectrum
                para_int_p, perp_int_p, freq, power_s_para, power_s_perp_1, power_s_perp_2 = FFT_Hann(mask)
                #save freq and power_s as their own df, named after 1- window start and 2- sc
                spectral_df = pd.DataFrame({'Freq':freq, 'Parallel Power':power_s_para, 'Perp 1 Power': power_s_perp_1, 'Perp 2 Power':power_s_perp_2})
                spec_fname = str(m)+'FS_C2.csv'
                fpath_spec = '/data/scratch/apx059/23_Years_Data/CSVs/GIPM_4mins/Fourier_Products/Fourier_Spectra_'+year_n+'/'+spec_fname
                spectral_df.to_csv(fpath_spec)
                
                cl_min_list.append(Cluster_min_rat)
                cl_mean_list.append(Cluster_mean_rat)
                cl_max_list.append(Cluster_max_rat)
                cl_median_list.append(Cluster_median_rat)
                cl_std_list.append(Cluster_stf_rat)
                cone_angle_list.append(omni_ave_c_a)
                para_i_p_list.append(para_int_p)
                perp_i_p_list.append(perp_int_p)
                ma_list.append(omni_ave_m_a)
                times.append(m)

        if cl_min_list:
            B_val_df = pd.DataFrame({'datetime': times,'B min': cl_min_list, 'B mean': cl_mean_list, 'B max': cl_max_list,'B median': cl_median_list,'B standard deviation':cl_std_list, 'cone angle': cone_angle_list, 'M_A': ma_list, 'ULF Parallel Power':para_i_p_list, 'ULF Perpendicular Power': perp_i_p_list})
            B_val_df = B_val_df.set_index('datetime')
            new_cl_df = j.join([B_val_df])
            list_expanded_dfs_2.append(new_cl_df)

        
    #for df in Cluster_GIPM_locs_list

    for df in list_expanded_dfs_2:
        if df.size > 0:
            firstwin = df.index[0]
            firstwin = str(firstwin)
            fpath = CSV_path + firstwin + 'C2.csv'
            df.to_csv(fpath)
            
    #for df in omni_ave list:

    for om_df in om_ave_list_2:
        if om_df.size > 0:
            firstwin = om_df.index[0]
            firstwin = str(firstwin)
            fpath = OMNI_CSV_path + firstwin + 'OMNI_C2.csv'
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
        cl_median_list = []
        cl_std_list = []
        cone_angle_list = []
        para_i_p_list = []
        perp_i_p_list = []
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
                Cluster_median = statistics.median(Cluster_list)
                Cluster_mean = sum(Cluster_list)/len(Cluster_list)
                Cluster_stf = statistics.stdev(Cluster_list)
                Cluster_stf_rat = Cluster_stf/omni_ave_B
                Cluster_min_rat = Cluster_min/omni_ave_B
                Cluster_mean_rat = Cluster_mean/omni_ave_B
                Cluster_max_rat = Cluster_max/omni_ave_B
                Cluster_median_rat = Cluster_median/omni_ave_B
                
                #4 minute power spectrum
                para_int_p, perp_int_p, freq, power_s_para, power_s_perp_1, power_s_perp_2 = FFT_Hann(mask)
                #save freq and power_s as their own df, named after 1- window start and 2- sc
                spectral_df = pd.DataFrame({'Freq':freq, 'Parallel Power':power_s_para, 'Perp 1 Power': power_s_perp_1, 'Perp 2 Power':power_s_perp_2})
                spec_fname = str(m)+'FS_C3.csv'
                fpath_spec = '/data/scratch/apx059/23_Years_Data/CSVs/GIPM_4mins/Fourier_Products/Fourier_Spectra_'+year_n+'/'+spec_fname
                spectral_df.to_csv(fpath_spec)
                
                
                cl_min_list.append(Cluster_min_rat)
                cl_mean_list.append(Cluster_mean_rat)
                cl_max_list.append(Cluster_max_rat)
                cl_median_list.append(Cluster_median_rat)
                cl_std_list.append(Cluster_stf_rat)
                cone_angle_list.append(omni_ave_c_a)
                para_i_p_list.append(para_int_p)
                perp_i_p_list.append(perp_int_p)
                ma_list.append(omni_ave_m_a)
                times.append(m)

        if cl_min_list:
            B_val_df = pd.DataFrame({'datetime': times,'B min': cl_min_list, 'B mean': cl_mean_list, 'B max': cl_max_list,'B median': cl_median_list, 'B standard deviation': cl_std_list, 'cone angle': cone_angle_list, 'M_A': ma_list, 'ULF Parallel Power':para_i_p_list, 'ULF Perpendicular Power': perp_i_p_list})
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
            fpath = OMNI_CSV_path + firstwin + 'OMNI_C3.csv'
            om_df.to_csv(fpath)
            
    ######C4
    
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
    GIPM_X_vec_list_4 = []
    GIPM_Y_vec_list_4 = []
    GIPM_Z_vec_list_4 = []
    FAC_coeff_list_4 = []

    for om_df in om_ave_list_4:
        GIPM_X_Vecs, GIPM_Y_Vecs, GIPM_Z_Vecs, FAC_coeffs = gipm_transform_coeffs(om_df)
        GIPM_X_vec_list_4.append(GIPM_X_Vecs)
        GIPM_Y_vec_list_4.append(GIPM_Y_Vecs)
        GIPM_Z_vec_list_4.append(GIPM_Z_Vecs)
        FAC_coeff_list_4.append(FAC_coeffs)

    Cluster_GIPM_locs_list_4 = []

    for p,q,i,j,k,m in zip(f_winds_all_4, df_list_c4, GIPM_X_vec_list_4, GIPM_Y_vec_list_4, GIPM_Z_vec_list_4, FAC_coeff_list_4):
        Cluster_dt_loc = gipm_loc_transform(p,q,i,j,k,m)
        Cluster_GIPM_locs_list_4.append(Cluster_dt_loc)
        
    gipm_df_c4 = []

    for df_c4 in Cluster_GIPM_locs_list_4:
        df_c4['datetime'] = pd.to_datetime(df_c4['datetime'])
        df_c4 = df_c4.set_index('datetime')
        gipm_df_c4.append(df_c4)

    #this bit is for the new dfs 

    list_expanded_dfs_4 = []

    time_window = dt.timedelta(seconds=240)

    for i,j,k in zip(df_list_c4, gipm_df_c4, om_ave_list_4):

        cl_min_list = []
        cl_max_list = []
        cl_mean_list = []
        cl_median_list = []
        cl_std_list = []
        cone_angle_list = []
        para_i_p_list = []
        perp_i_p_list = []
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
                Cluster_median = statistics.median(Cluster_list)
                Cluster_mean = sum(Cluster_list)/len(Cluster_list)
                Cluster_stf = statistics.stdev(Cluster_list)
                Cluster_stf_rat = Cluster_stf/omni_ave_B
                Cluster_min_rat = Cluster_min/omni_ave_B
                Cluster_mean_rat = Cluster_mean/omni_ave_B
                Cluster_max_rat = Cluster_max/omni_ave_B
                Cluster_median_rat = Cluster_median/omni_ave_B
                
                #4 minute power spectrum
                para_int_p, perp_int_p, freq, power_s_para, power_s_perp_1, power_s_perp_2 = FFT_Hann(mask)
                #save freq and power_s as their own df, named after 1- window start and 2- sc
                spectral_df = pd.DataFrame({'Freq':freq, 'Parallel Power':power_s_para, 'Perp 1 Power': power_s_perp_1, 'Perp 2 Power':power_s_perp_2})
                spec_fname = str(m)+'_FS_C4.csv'
                fpath_spec = '/data/scratch/apx059/23_Years_Data/CSVs/GIPM_4mins/Fourier_Products/Fourier_Spectra_'+year_n+'/'+spec_fname
                spectral_df.to_csv(fpath_spec)
                
                cl_min_list.append(Cluster_min_rat)
                cl_mean_list.append(Cluster_mean_rat)
                cl_max_list.append(Cluster_max_rat)
                cl_median_list.append(Cluster_median_rat)
                cl_std_list.append(Cluster_stf_rat)
                cone_angle_list.append(omni_ave_c_a)
                para_i_p_list.append(para_int_p)
                perp_i_p_list.append(perp_int_p)
                ma_list.append(omni_ave_m_a)
                times.append(m)
                
        if cl_min_list:
            B_val_df = pd.DataFrame({'datetime': times,'B min': cl_min_list, 'B mean': cl_mean_list, 'B max': cl_max_list,'B median': cl_median_list, 'B standard deviation':cl_std_list, 'cone angle': cone_angle_list, 'M_A': ma_list, 'ULF Parallel Power':para_i_p_list, 'ULF Perpendicular Power': perp_i_p_list})
            B_val_df = B_val_df.set_index('datetime')
            new_cl_df = j.join([B_val_df])
            list_expanded_dfs_4.append(new_cl_df)

        
    #for df in Cluster_GIPM_locs_list
    CSV_path = '/data/scratch/apx059/23_Years_Data/CSVs/GIPM_4mins/'+year_n+'_Integrated_CSVs/'
    OMNI_CSV_path = '/data/scratch/apx059/23_Years_Data/CSVs/GIPM_4mins/'+year_n+'_OMNI_CSVs/'

    for df in list_expanded_dfs_4:
        if df.size > 0:
            firstwin = df.index[0]
            firstwin = str(firstwin)
            fpath = CSV_path + firstwin + 'C4.csv'
            df.to_csv(fpath)
            
    #for df in omni_ave list:

    for om_df in om_ave_list_4:
        if om_df.size > 0:
            firstwin = om_df.index[0]
            firstwin = str(firstwin)
            fpath = OMNI_CSV_path + firstwin + 'OMNI_C4.csv'
            om_df.to_csv(fpath)
    
    #########################

cdfconv_gipm(year_n)


