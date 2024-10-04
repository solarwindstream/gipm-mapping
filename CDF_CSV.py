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
from gipm_locs_quick import gipm_locs_quick

year_n = '2002'
cdf_list_path = '/data/scratch/apx059/23_Years_Data/CDFs/CDFs_'
list_of_cdfs = cdf_list_path + year_n + '.csv'

def cdfconv_gipm(path, year):
    
    #batch start and end
    batch_start = batch*100
    batch_end = (batch+1)*100
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
            df_list_c1.append(df_c1)
        if 'C2' in i:
            fpath = year_path + i
            df_c2 = C2_cdf_conv(fpath)
            df_list_c2.append(df_c2)
        if 'C3' in i:   
            fpath = year_path + i
            df_c3 = C3_cdf_conv(fpath)
            df_list_c3.append(df_c3)   
        if 'C4' in i:   
            fpath = year_path + i
            df_c4 = C4_cdf_conv(fpath)
            df_list_c4.append(df_c4)
            
    #first append dataframes that *aren't* empty to new list (to avoid errors)
    df_full_c1 = []
    df_full_c2 = []
    df_full_c3 = []
    df_full_c4 = []

    for df in df_list_c1:
        a = df.empty
        if not a:
            df_full_c1.append(df)

    for df in df_list_c2:
        a = df.empty
        if not a:
            df_full_c2.append(df)

    for df in df_list_c3:
        a = df.empty
        if not a:
            df_full_c3.append(df)

    for df in df_list_c4:
        a = df.empty
        if not a:
            df_full_c4.append(df) 


    #list all OMNI files, but only import relevant ones
    #load ONLY this year's omni df
    
    omni_csv_list_path = r'/data/home/apx059/gipm-mapping/listofomnis.csv'
    omni_files = pd.read_csv(omni_csv_list_path, header=None)
    omni_files = omni_files.rename(columns={0:'fnames'})
    list_omni_csvs = omni_files['fnames'].to_list()
    res = [i for i in list_omni_csvs if year_n in i]
    om = pd.read_csv(res[0])

    om['datetime'] = pd.to_datetime(om['datetime'])
    om = om.set_index('datetime')

    ##########################
    #GIPM conversion!
    ####PRESERVE SC LISTINGS
    ####SC1
    
    om_ave_list_1 = []

    for df in df_full_c1:
        datetime_list = df.index
        om_averages = omni_seg(om, datetime_list)
        om_ave_list_1.append(om_averages)

    for om_ave_df in om_ave_list_1:
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
        filename = CSV_path + 'Feb' + start_year + '_Feb' + end_year + 'df_' + ref + '.csv'
        np.savetxt(filename,array_both, fmt='%s,%f,%f,%f')
        startref = startref + 1

    #########################
    
csvconv(list_of_cdfs, year_n)