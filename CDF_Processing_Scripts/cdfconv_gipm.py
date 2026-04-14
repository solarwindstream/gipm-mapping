import datetime as dt
import pandas as pd
import numpy as np
import statistics

from window_det_4min import window_det
from omni_seg_centered import omni_seg, omni_seg_ap_ratio
from gipm_transform_coeffs import gipm_transform_coeffs_mean, gipm_transform_coeffs_median
from GIPM_loc_conv import gipm_loc_transform
from new_xyz import new_xyz
from FFT_Hann import FFT_Hann

#define the main conversion module!

def cdfconv_gipm(year, df_list, sc_name):
    
    #determine full windows lists.
    f_winds_all = []

    for df in df_list:
        f_winds = window_det(df)
        f_winds_all.append(f_winds)
    
    om_ave_list = []

    #calculate B,V,n,MA averages separately to Na/Np averages as they're saved separately, then join.
    for fw in f_winds_all:
        om_averages = omni_seg(om, fw)
        om_nanp_averages = omni_seg_ap_ratio(omni_nanp_df, fw)
        om_full_averages = om_averages.join(om_nanp_averages)
        om_ave_list.append(om_full_averages)

    #print('Omni averages found')

    #find GIPM rotation matrices and scaling coefficient for every Cluster location
    
    GIPM_X_mean_list = []
    GIPM_Y_mean_list = []
    GIPM_Z_mean_list = []
    FAC_coeff_mean_list = []
    
    GIPM_X_median_list = []
    GIPM_Y_median_list = []
    GIPM_Z_median_list = []
    FAC_coeff_median_list = []

    for om_df in om_ave_list:
        
        #mean
        GIPM_X_Vecs, GIPM_Y_Vecs, GIPM_Z_Vecs, FAC_coeffs = gipm_transform_coeffs_mean(om_df)
        GIPM_X_mean_list.append(GIPM_X_Vecs)
        GIPM_Y_mean_list.append(GIPM_Y_Vecs)
        GIPM_Z_mean_list.append(GIPM_Z_Vecs)
        FAC_coeff_mean_list.append(FAC_coeffs)
        
        #median
        GIPM_X_Vecs, GIPM_Y_Vecs, GIPM_Z_Vecs, FAC_coeffs = gipm_transform_coeffs_median(om_df)
        GIPM_X_median_list.append(GIPM_X_Vecs)
        GIPM_Y_median_list.append(GIPM_Y_Vecs)
        GIPM_Z_median_list.append(GIPM_Z_Vecs)
        FAC_coeff_median_list.append(FAC_coeffs)

    
    #now combine these basis vectors and co-efficients with Cluster GSE location data to transform into GIPM
    Cluster_GIPM_locs_list = []
    
    #mean AND median locations both produced by gipm_loc_transform
    for fw,df,GIPM_Xax_mean,GIPM_Yax_mean,GIPM_Zax_mean,FAC_mean,GIPM_Xax_median,GIPM_Yax_median,GIPM_Zax_median,FAC_median in zip(f_winds_all, df_list, GIPM_X_mean_list, GIPM_Y_mean_list, GIPM_Z_mean_list, FAC_coeff_mean_list, GIPM_X_median_list, GIPM_Y_median_list, GIPM_Z_median_list, FAC_coeff_median_list):
        Cluster_dt_loc_GIPM = gipm_loc_transform(fw,df,GIPM_Xax_mean,GIPM_Yax_mean,GIPM_Zax_mean,FAC_mean,GIPM_Xax_median,GIPM_Yax_median,GIPM_Zax_median,FAC_median)
        Cluster_GIPM_locs_list.append(Cluster_dt_loc_GIPM)
        
    #print('GIPM transform done')
        
    gipm_dfs = []

    for df in Cluster_GIPM_locs_list:
        df['datetime'] = pd.to_datetime(df['datetime'])
        df = df.set_index('datetime')
        gipm_dfs.append(df)

    #now find min, max and mean B, integrated power, and save with relevant OMNI average values.

    list_expanded_dfs = []

    time_window = dt.timedelta(seconds=240)

    for df,gipm_df,omni_ave_df in zip(df_list, gipm_dfs, om_ave_list):

        #empty lists set up for cluster variables
        cl_min_list, cl_mean_list, cl_max_list, cl_median_list, cl_std_list, para_i_p_list, perp_i_p_list, para_norm_i_p_list, perp_norm_i_p_list, peak_comp_freq_list, peak_trans_freq_list, ellipticity_list, times = ([] for i in range(13))
        #empty lists set up for core upstream variables
        cone_angle_mean_list, cone_angle_med_list, ma_mean_list, ma_med_list, omni_mean_B_list, omni_med_B_list, omni_mean_V_list, omni_med_V_list, omni_mean_Np_list, omni_med_Np_list, omni_ave_max_IMF_dev_list, omni_sc_dist_mean_list, omni_sc_dist_median_list = ([] for i in range(13))
        #empty lists for additional upstream variables (mean and median B and V vecs, NaNp ratios)
        omni_mean_Bx_list, omni_mean_By_list, omni_mean_Bz_list, omni_mean_Vx_list, omni_mean_Vy_list, omni_mean_Vz_list, omni_med_Bx_list, omni_med_By_list, omni_med_Bz_list, omni_med_Vx_list, omni_med_Vy_list, omni_med_Vz_list, omni_mean_NaNp_list, omni_med_NaNp_list = ([] for i in range(14))
        
        #drop any intervals without GIPM location values!
        gipm_df.dropna(subset=['GIPM X (OMNI mean)'], inplace=True)
        
        if gipm_df.size > 0:
        
            for win_start in gipm_df.index:
                start_time = win_start
                end_time = win_start + time_window
                mask = df.loc[(df.index >= start_time) & (df.index < end_time)]
                Cluster_list = mask['B_mag'].tolist()
                
                omni_ave_B = omni_ave_df.loc[win_start, 'B_mag (mean)']
                omni_med_B = omni_ave_df.loc[win_start, 'B_mag (median)']
                omni_ave_Bx = omni_ave_df.loc[win_start, 'B_X_gse (mean)']
                omni_ave_By = omni_ave_df.loc[win_start, 'B_Y_gse (mean)']
                omni_ave_Bz = omni_ave_df.loc[win_start, 'B_Z_gse (mean)']
                omni_med_Bx = omni_ave_df.loc[win_start, 'B_X_gse (median)']
                omni_med_By = omni_ave_df.loc[win_start, 'B_Y_gse (median)']
                omni_med_Bz = omni_ave_df.loc[win_start, 'B_Z_gse (median)']
                omni_ave_V = omni_ave_df.loc[win_start, 'V_gse (mean)']
                omni_med_V = omni_ave_df.loc[win_start, 'V_gse (median)']
                omni_ave_Vx = omni_ave_df.loc[win_start, 'V_X_gse (mean)']
                omni_ave_Vy = omni_ave_df.loc[win_start, 'V_Y_gse (mean)']
                omni_ave_Vz = omni_ave_df.loc[win_start, 'V_Z_gse (mean)']
                omni_med_Vx = omni_ave_df.loc[win_start, 'V_X_gse (median)']
                omni_med_Vy = omni_ave_df.loc[win_start, 'V_Y_gse (median)']
                omni_med_Vz = omni_ave_df.loc[win_start, 'V_Z_gse (median)']
                omni_ave_Np = omni_ave_df.loc[win_start, 'Np (mean)']
                omni_med_Np = omni_ave_df.loc[win_start, 'Np (median)']
                omni_ave_NaNp = omni_ave_df.loc[win_start, 'Na_Np (mean)']
                omni_med_NaNp = omni_ave_df.loc[win_start, 'Na_Np (median)']
                omni_ave_c_a = omni_ave_df.loc[win_start, 'cone angle (mean)']
                omni_med_c_a = omni_ave_df.loc[win_start, 'cone angle (median)']
                omni_ave_m_a = omni_ave_df.loc[win_start, 'M_A (mean)']
                omni_med_m_a = omni_ave_df.loc[win_start, 'M_A (median)']
                omni_ave_max_IMF_dev = omni_ave_df.loc[win_start, 'max IMF angle deviation']
                omni_sc_dist_mean = omni_ave_df.loc[win_start, 'Distance from X line (mean)']
                omni_sc_dist_median = omni_ave_df.loc[win_start, 'Distance from X line (median)']
                
                
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
                    
                    #calculate power spectrum and integrated power for this window
                    para_int_p, perp_int_p, ellipticity, freq, power_s_para, power_s_perp_1, power_s_perp_2 = FFT_Hann(mask)
                    
                    #save the power spectra as their own df, filename including window start time & sc number
                    spectral_df = pd.DataFrame({'Frequency':freq, 'Compressive Power':power_s_para, 'Transverse Power 1': power_s_perp_1, 'Transverse Power 2':power_s_perp_2})
                    spectral_df['Total Transverse Power'] = spectral_df[ 'Transverse Power 1'] + spectral_df[ 'Transverse Power 2']
                    spec_fname = 'FS_' + str(win_start) + '_' + sc_name + '.csv'
                    fpath_spec = '/Users/roseatkinson/Documents/Debug Output/Spectra_'+ year + '/' + spec_fname
                    spectral_df.to_csv(fpath_spec)

                    #find the index of peak compressive power, zero removed
                    spectral_df_no_zero = spectral_df.iloc[1:]
                    max_comp_pwr_idx = spectral_df_no_zero[['Compressive Power']].idxmax()
                    max_freq_comp = spectral_df_no_zero.loc[max_comp_pwr_idx, 'Frequency'].values[0]
                    max_trans_pwr_idx = spectral_df_no_zero[['Total Transverse Power']].idxmax()
                    max_freq_trans = spectral_df_no_zero.loc[max_trans_pwr_idx, 'Frequency'].values[0]
                    
                    #normalise power according to B^2:
                    para_int_p_norm = para_int_p/(Cluster_mean**2)
                    perp_int_p_norm = perp_int_p/(Cluster_mean**2)

                    cl_min_list.append(Cluster_min_rat)
                    cl_mean_list.append(Cluster_mean_rat)
                    cl_max_list.append(Cluster_max_rat)
                    cl_median_list.append(Cluster_median_rat)
                    cl_std_list.append(Cluster_stf_rat)

                    #add new OMNI values
                    cone_angle_mean_list.append(omni_ave_c_a)
                    cone_angle_med_list.append(omni_med_c_a)
                    ma_mean_list.append(omni_ave_m_a)
                    ma_med_list.append(omni_med_m_a)
                    omni_mean_B_list.append(omni_ave_B)
                    omni_med_B_list.append(omni_med_B)
                    omni_mean_V_list.append(omni_ave_V)
                    omni_med_V_list.append(omni_med_V)
                    omni_mean_Np_list.append(omni_ave_Np)
                    omni_med_Np_list.append(omni_med_Np)

                    omni_mean_Bx_list.append(omni_ave_Bx)
                    omni_mean_By_list.append(omni_ave_By)
                    omni_mean_Bz_list.append(omni_ave_Bz)
                    omni_mean_Vx_list.append(omni_ave_Vx)
                    omni_mean_Vy_list.append(omni_ave_Vy)
                    omni_mean_Vz_list.append(omni_ave_Vz)
                    omni_med_Bx_list.append(omni_med_Bx)
                    omni_med_By_list.append(omni_med_By)
                    omni_med_Bz_list.append(omni_med_Bz)
                    omni_med_Vx_list.append(omni_med_Vx)
                    omni_med_Vy_list.append(omni_med_Vy)
                    omni_med_Vz_list.append(omni_med_Vz)

                    omni_mean_NaNp_list.append(omni_ave_NaNp)
                    omni_med_NaNp_list.append(omni_med_NaNp)
                    
                    omni_ave_max_IMF_dev_list.append(omni_ave_max_IMF_dev)
                    omni_sc_dist_mean_list.append(omni_sc_dist_mean)
                    omni_sc_dist_median_list.append(omni_sc_dist_median)
                    
                    #add spectral values 
                    
                    para_i_p_list.append(para_int_p)
                    perp_i_p_list.append(perp_int_p)
    
                    para_norm_i_p_list.append(para_int_p_norm)
                    perp_norm_i_p_list.append(perp_int_p_norm)

                    peak_comp_freq_list.append(max_freq_comp)
                    peak_trans_freq_list.append(max_freq_trans)
                    ellipticity_list.append(ellipticity)
                    
                    times.append(win_start)
    
                else:
                    print('error: no cluster list')
            
            if cl_min_list:
                core_val_df = pd.DataFrame({'datetime': times,'B min/Bomni': cl_min_list, 'B mean/Bomni': cl_mean_list, 'B max/Bomni': cl_max_list, 'B median/Bomni': cl_median_list, 'B standard deviation/Bomni':cl_std_list, 'cone angle (mean)': cone_angle_mean_list, 'M_A (mean)': ma_mean_list, 'IMF B (mean)': omni_mean_B_list, 'SW V (mean)': omni_mean_V_list, 'SW Np (mean)': omni_mean_Np_list, 'SW Na/Np (mean)': omni_mean_NaNp_list, 'OMNI Dist from X line (mean)': omni_sc_dist_mean_list, 'Max IMF Deviation': omni_ave_max_IMF_dev_list, 'ULF Band Compressive Power': para_i_p_list, 'ULF Band Transverse Power': perp_i_p_list, 'ULF Band Normalised Compressive Power': para_norm_i_p_list, 'ULF Band Normalised Transverse Power': perp_norm_i_p_list, 'Peak Compressive Frequency': peak_comp_freq_list, 'Peak Transverse Frequency': peak_trans_freq_list,'Ratio of Perpendicular Power': ellipticity_list})
                core_val_df = core_val_df.set_index('datetime')
                secondary_val_df = pd.DataFrame({'datetime': times, 'cone angle (median)': cone_angle_med_list, 'M_A (median)': ma_med_list, 'IMF B (median)': omni_med_B_list, 'SW V (median)': omni_med_V_list, 'SW Np (median)': omni_med_Np_list, 'SW Na/Np (median)': omni_med_NaNp_list, 'OMNI Dist from X line (median)': omni_sc_dist_median_list, 'IMF Bx (mean)': omni_mean_Bx_list, 'IMF By (mean)': omni_mean_By_list, 'IMF Bz (mean)': omni_mean_Bz_list, 'SW Vx (mean)': omni_mean_Vx_list, 'SW Vy (mean)': omni_mean_Vy_list, 'SW Vz (mean)': omni_mean_Vz_list, 'IMF Bx (median)': omni_med_Bx_list, 'IMF By (median)': omni_med_By_list, 'IMF Bz (median)': omni_med_Bz_list, 'SW Vx (median)': omni_med_Vx_list, 'SW Vy (median)': omni_med_Vy_list, 'SW Vz (median)': omni_med_Vz_list})
                secondary_val_df = secondary_val_df.set_index('datetime')
                new_cl_df = gipm_df.join([core_val_df, secondary_val_df])
                list_expanded_dfs.append(new_cl_df)
            else:
                print('error: no cluster min list')
            
    #for df in Cluster_GIPM_locs_list
    CSV_path = '/Users/roseatkinson/Documents/Debug Output/'

    for df in list_expanded_dfs:
        firstwin = df.index[0]
        firstwin = str(firstwin)
        fpath = CSV_path + firstwin + '_' + sc_name + '.csv'
        df.to_csv(fpath)
