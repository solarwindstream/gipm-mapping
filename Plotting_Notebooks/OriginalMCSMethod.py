import pandas as pd
import numpy as np
import datetime as dt
import glob

import matplotlib 
from matplotlib import pyplot as plt
from spectra_modules import spectra_plot, TS_plot, spectra_overlay_plot, TS_and_FS, split_band_power_int, spectra_overlay_plot_comparison
from gipm_plotting_modules import single_map, location_mapping
from heatmap_modules import compute_hists2d
import scipy
##load integrated Cluster CSVs

cl_file_list_C1 = []
cl_file_list_C2 = []
cl_file_list_C3 = []
cl_file_list_C4 = []

path = "/Users/roseatkinson/Documents/Cluster_Int_CSVs/**"

for path in glob.glob(path, recursive=True):
    if '.csv' in path:
        if 'C1' in path:
            cl_file_list_C1.append(path)
        if 'C2' in path:
            cl_file_list_C2.append(path)
        if 'C3' in path:
            cl_file_list_C3.append(path)
        if 'C4' in path:
            cl_file_list_C4.append(path)
    
cl_dfs = []
sc_list = ['C1', 'C2', 'C3', 'C4']
all_lists = [cl_file_list_C1, cl_file_list_C2, cl_file_list_C3, cl_file_list_C4]

for sc, file_list in zip(sc_list, all_lists):
    for file in file_list:
        df = pd.read_csv(file,encoding='utf-8')
        df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
        df.set_index('datetime', inplace = True)
        df['sc']=sc
        cl_dfs.append(df)

cl_power_all = pd.concat(cl_dfs)

print("dfs loaded in")

#Filtering: remove times when OMNI source sc was too far from X-line, & current sheets, & filter down to -5<Zgipm<5 "plane"
cl_filtered = cl_power_all.loc[(cl_power_all['OMNI Dist from X line (mean)'] < 70) & (cl_power_all['Max IMF Deviation'] < 60)]
cl_filtered_lowZ = cl_filtered.loc[(cl_filtered['GIPM Z (OMNI mean)'] < 5) & (cl_filtered['GIPM Z (OMNI mean)'] > -5)]

##load extreme time series Cluster CSVs

cl_ex_dfs_C1 = []
cl_ex_dfs_C2 = []
cl_ex_dfs_C3 = []
cl_ex_dfs_C4 = []

path = "/Users/roseatkinson/Documents/Cluster_Time_Series/**"

for path in glob.glob(path, recursive=True):
    if '.csv' in path:
        df = pd.read_csv(path,encoding='utf-8')
        df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
        df.set_index('datetime', inplace = True)
        if 'C1' in path:
            cl_ex_dfs_C1.append(df)
        if 'C2' in path:
            cl_ex_dfs_C2.append(df)
        if 'C3' in path:
            cl_ex_dfs_C3.append(df)
        if 'C4' in path:
            cl_ex_dfs_C4.append(df)

cl_ts_1 = pd.concat(cl_ex_dfs_C1)
cl_ts_2 = pd.concat(cl_ex_dfs_C2)
cl_ts_3 = pd.concat(cl_ex_dfs_C3)
cl_ts_4 = pd.concat(cl_ex_dfs_C4)

cl_ts_1.sort_index(inplace=True)
cl_ts_2.sort_index(inplace=True)
cl_ts_3.sort_index(inplace=True)
cl_ts_4.sort_index(inplace=True)

#filtered for 3<B<6
cl_filtered_lowZ = cl_filtered_lowZ.loc[(cl_filtered_lowZ['IMF B (mean)'] > 3) & (cl_filtered_lowZ['IMF B (mean)'] < 6)]
cl_filtered_40_60 = cl_filtered_lowZ.loc[(cl_filtered_lowZ['cone angle (mean)'] >= 40) & (cl_filtered_lowZ['cone angle (mean)'] < 60)]
cl_filtered_20_40 = cl_filtered_lowZ.loc[(cl_filtered_lowZ['cone angle (mean)'] >= 20) & (cl_filtered_lowZ['cone angle (mean)'] < 40)]
cl_filtered_0_20 = cl_filtered_lowZ.loc[(cl_filtered_lowZ['cone angle (mean)'] >= 0) & (cl_filtered_lowZ['cone angle (mean)'] < 20)]

heatmaps_dict_new = {}
heatmaps_dict_new['0-20°'] = {}
heatmaps_dict_new['20-40°'] = {}
heatmaps_dict_new['40-60°'] = {}

heatmaps_dict_new['0-20°'] = compute_hists2d(cl_filtered_0_20, 'mean', 'Normalised Transverse Power', obs_min=50)
heatmaps_dict_new['20-40°'] = compute_hists2d(cl_filtered_20_40, 'mean', 'Normalised Transverse Power', obs_min=50)
heatmaps_dict_new['40-60°'] = compute_hists2d(cl_filtered_40_60, 'mean', 'Normalised Transverse Power', obs_min=50)

#note that ion and tangent lines are for 30-45 case
#single_map(heatmaps_dict_new['0-20°']['Normalised Transverse Power'], "0–20°",'Power')
#single_map(heatmaps_dict_new['20-40°']['Normalised Transverse Power'], "20–40°",'Power')
#single_map(heatmaps_dict_new['40-60°']['Normalised Transverse Power'], "40–60°",'Power')

#filter for likely ULF
cl_filtered_0_20_ULF = cl_filtered_0_20.loc[(cl_filtered_0_20['ULF Band Normalised Transverse Power'] > 0.05)]
cl_filtered_20_40_ULF = cl_filtered_20_40.loc[(cl_filtered_20_40['ULF Band Normalised Transverse Power'] > 0.05)]
cl_filtered_40_60_ULF = cl_filtered_40_60.loc[(cl_filtered_40_60['ULF Band Normalised Transverse Power'] > 0.05)]

#now iterate over x, y foreshock strips:
#key is start of x pair

y_lims_0_20 = {
    '9':[-16,-12], 
    '10':[-16,-11], 
    '11':[-15,-8], 
    '12':[-14,-6], 
    '13':[-13,8],
    '14':[-12,8],
    '15':[-11,3]
}

y_lims_20_40 = {
    '5':[-18,-16], 
    '6':[-18,-16], 
    '7':[-17,-15],
    '8':[-16,-13], 
    '9':[-16,-10], 
    '10':[-16,-9], 
    '11':[-14,-6], 
    '12':[-13,-6], 
    '13':[-13,1],
    '14':[-13,0],
    '15':[-12,0],
    '16':[-10,0]
}

y_lims_40_60 = {
    '6':[-18,-16], 
    '7':[-17,-14], 
    '8':[-17,-12], 
    '9':[-17,-12], 
    '10':[-17,-12],
    '11':[-16,-9],
    '12':[-14,-7],
}

fs_df_list_0_20 = []
fs_df_list_20_40 = []
fs_df_list_40_60 = []

for x in range(9,16):
    x_end = x + 1
    y_bounds = y_lims_0_20[str(x)]
    df = cl_filtered_0_20_ULF.loc[(cl_filtered_0_20_ULF['GIPM X (OMNI mean)'] >= x) & (cl_filtered_0_20_ULF['GIPM X (OMNI mean)'] < x_end) & (cl_filtered_0_20_ULF['GIPM Y (OMNI mean)'] >= y_bounds[0]) & (cl_filtered_0_20_ULF['GIPM Y (OMNI mean)'] < y_bounds[1])]
    fs_df_list_0_20.append(df)

for x in range(5,17):
    x_end = x + 1
    y_bounds = y_lims_20_40[str(x)]
    df = cl_filtered_20_40_ULF.loc[(cl_filtered_20_40_ULF['GIPM X (OMNI mean)'] >= x) & (cl_filtered_20_40_ULF['GIPM X (OMNI mean)'] < x_end) & (cl_filtered_20_40_ULF['GIPM Y (OMNI mean)'] >= y_bounds[0]) & (cl_filtered_20_40_ULF['GIPM Y (OMNI mean)'] < y_bounds[1])]
    fs_df_list_20_40.append(df)

for x in range(6,13):
    x_end = x + 1
    y_bounds = y_lims_40_60[str(x)]
    df = cl_filtered_40_60_ULF.loc[(cl_filtered_40_60_ULF['GIPM X (OMNI mean)'] >= x) & (cl_filtered_40_60_ULF['GIPM X (OMNI mean)'] < x_end) & (cl_filtered_40_60_ULF['GIPM Y (OMNI mean)'] >= y_bounds[0]) & (cl_filtered_40_60_ULF['GIPM Y (OMNI mean)'] < y_bounds[1])]
    fs_df_list_40_60.append(df)

fs_dfs_0_20 = pd.concat(fs_df_list_0_20)
fs_dfs_20_40 = pd.concat(fs_df_list_20_40)
fs_dfs_40_60 = pd.concat(fs_df_list_40_60)

cl_fs_highNaNp_0_20 = fs_dfs_0_20.loc[(fs_dfs_0_20['SW Na/Np (mean)'] > 0.08)]
cl_fs_lowNaNp_0_20 = fs_dfs_0_20.loc[(fs_dfs_0_20['SW Na/Np (mean)'] < 0.005)]
cl_fs_highNaNp_20_40 = fs_dfs_20_40.loc[(fs_dfs_20_40['SW Na/Np (mean)'] > 0.08)]
cl_fs_lowNaNp_20_40 = fs_dfs_20_40.loc[(fs_dfs_20_40['SW Na/Np (mean)'] < 0.005)]
cl_fs_highNaNp_40_60 = fs_dfs_40_60.loc[(fs_dfs_40_60['SW Na/Np (mean)'] > 0.08)]
cl_fs_lowNaNp_40_60 = fs_dfs_40_60.loc[(fs_dfs_40_60['SW Na/Np (mean)'] < 0.005)]

#low Nanp, 0-20deg. 
path = '/Users/roseatkinson/Documents/Extreme_Spectra/'

fs_df_list_low_20 = []
win_starts_low_20 = []
sc_no_low_20 = []
labels_low_20 = []
omni_means_low_20 = []
x_means_low_20 = []
y_means_low_20 = []
error_files_low_20 = []
takahashi_low_20 = []
mean_B_cluster_low_20 = []
mag_comp_low_20 = []
v_sw_means_low_20 = []
cone_angle_low_20 = []
proton_density_low_20 = []

for datetime, row in cl_fs_lowNaNp_0_20.iterrows():
    int_start = datetime.strftime("%Y-%m-%d %H:%M:%S")
    sc_no = row['sc']
    filename = 'FS_'+int_start+'_'+sc_no+'.csv'
    filepath = path+filename
    try:
        f_df = pd.read_csv(filepath,encoding='utf-8')
    except:
        error_files_low_20.append(filename)
    else:
        fs_df_list_low_20.append(f_df)
        win_starts_low_20.append(datetime)
        alpha_comp = row['SW Na/Np (mean)']
        alpha_comp = round(alpha_comp, 3)
        label = r'$N\alpha /Np$=' + str(alpha_comp)
        labels_low_20.append(label)
        sc_no_low_20.append(sc_no)
        omni_mean = row['IMF B (mean)']
        omni_means_low_20.append(round(omni_mean,1))
        x_gipm_mean = row['GIPM X (OMNI mean)']
        y_gipm_mean = row['GIPM Y (OMNI mean)']
        tak_freq = row['Takahashi Frequency, Hz']
        takahashi_low_20.append(tak_freq)
        x_means_low_20.append(x_gipm_mean)
        y_means_low_20.append(y_gipm_mean)
        bmean_bomni = row['B mean/Bomni']
        mag_comp_low_20.append(bmean_bomni)
        cluster_b_mean = bmean_bomni*omni_mean
        mean_B_cluster_low_20.append(cluster_b_mean)
        v_sw_means_low_20.append(round(row['SW V (mean)']))
        cone_angle_low_20.append(round(row['cone angle (mean)']))
        proton_density_low_20.append(round(row['SW Np (mean)'], 1)) 
        
ts_dfs_low_20 = []

for win_start, sc in zip(win_starts_low_20, sc_no_low_20):
    win_end = win_start + dt.timedelta(minutes=4)
    sc_dict = {'C1': cl_ts_1, 'C2': cl_ts_2, 'C3': cl_ts_3, 'C4': cl_ts_4}
    ts_df = sc_dict[sc]
    ts_df_filt = ts_df.loc[((ts_df.index >= win_start) & (ts_df.index < win_end))]
    ts_dfs_low_20.append(ts_df_filt)

fs_dfs_normed_low_20 = []
int_perp_Pc4_low_20 = []
int_perp_Pc3_low_20 = []
int_para_Pc4_low_20 = []
int_para_Pc3_low_20 = []

for df, B in zip(fs_df_list_low_20, mean_B_cluster_low_20):
    new_df = df.copy()
    new_df['Total Transverse Power'] = new_df['Total Transverse Power']/(B**2)
    new_df['Compressive Power'] = new_df['Compressive Power']/(B**2)
    fs_dfs_normed_low_20.append(new_df)
    P_perp_low, P_perp_high, P_para_low, P_para_high = split_band_power_int(new_df)
    int_perp_Pc4_low_20.append(P_perp_low)
    int_perp_Pc3_low_20.append(P_perp_high)
    int_para_Pc4_low_20.append(P_para_low)
    int_para_Pc3_low_20.append(P_para_high)
    
#find time series
#from concatenated cluster dfs!
ratio_Pc4_Pc3_perp_low_20 = np.array(int_perp_Pc4_low_20)/np.array(int_perp_Pc3_low_20)
ratio_Pc4_Pc3_para_low_20 = np.array(int_para_Pc4_low_20)/np.array(int_para_Pc3_low_20)

# TS_and_FS(ts_dfs_low_20, fs_df_list_low_20, labels_low_20, tak_freq=takahashi_low_20)
#spectra_overlay_plot(fs_dfs_normed_low_20, 'low', '0-20°')

x_locs_low_0_20 = np.array(x_means_low_20)
y_locs_low_0_20 = np.array(y_means_low_20)

xy_pairs_low_0_20 = np.stack((x_locs_low_0_20,y_locs_low_0_20), axis=-1)

#location_mapping(xy_pairs_low_0_20, "0–20°", filename='low_0_20_cases')

#0-20
#relevant low Nanp: 0,1 in list.
#relevant high Nanp: 0-3 in list

#new set of lists

low_0_20_locations_A = []
low_0_20_ts_dfs_A = []
low_0_20_fs_dfs_normed_A = []
low_0_20_labels_A = []
low_0_20_tak_A = []
low_0_20_ratioPc4_3_perp_A = []
low_0_20_ratioPc4_3_para_A = []
low_0_20_V_SW_A = []
low_0_20_cone_angle_A = []
low_0_20_n_P_A = []
low_0_20_B_IMF_A = []

for ind in [0,1]:
    low_0_20_locations_A.append(xy_pairs_low_0_20[ind])
    low_0_20_ts_dfs_A.append(ts_dfs_low_20[ind])
    low_0_20_labels_A.append(labels_low_20[ind])
    low_0_20_tak_A.append(takahashi_low_20[ind])
    low_0_20_fs_dfs_normed_A.append(fs_dfs_normed_low_20[ind])
    low_0_20_ratioPc4_3_perp_A.append(ratio_Pc4_Pc3_perp_low_20[ind])
    low_0_20_ratioPc4_3_para_A.append(ratio_Pc4_Pc3_para_low_20[ind])
    low_0_20_V_SW_A.append(v_sw_means_low_20[ind])
    low_0_20_cone_angle_A.append(cone_angle_low_20[ind])
    low_0_20_n_P_A.append(proton_density_low_20[ind])
    low_0_20_B_IMF_A.append(omni_means_low_20[ind])

high_0_20_locations_A = []
high_0_20_ts_dfs_A = []
high_0_20_fs_dfs_normed_A = []
high_0_20_labels_A = []
high_0_20_tak_A = []
high_0_20_ratioPc4_3_perp_A = []
high_0_20_ratioPc4_3_para_A = []
high_0_20_V_SW_A = []
high_0_20_cone_angle_A = []
high_0_20_n_P_A = []
high_0_20_B_IMF_A = []

for ind in [0,1,2,3]:
    high_0_20_locations_A.append(xy_pairs_high_0_20[ind])
    high_0_20_ts_dfs_A.append(ts_dfs_high_20[ind])
    high_0_20_labels_A.append(labels_high_20[ind])
    high_0_20_tak_A.append(takahashi_high_20[ind])
    high_0_20_fs_dfs_normed_A.append(fs_dfs_normed_high_20[ind])
    high_0_20_ratioPc4_3_perp_A.append(ratio_Pc4_Pc3_perp_high_20[ind])
    high_0_20_ratioPc4_3_para_A.append(ratio_Pc4_Pc3_para_high_20[ind])
    high_0_20_V_SW_A.append(v_sw_means_high_20[ind])
    high_0_20_cone_angle_A.append(cone_angle_high_20[ind])
    high_0_20_n_P_A.append(proton_density_high_20[ind])
    high_0_20_B_IMF_A.append(omni_means_high_20[ind])

location_mapping(low_0_20_locations_A, "0–20°",  alt_locs=high_0_20_locations_A, filename='low_high_0_20_cases_A')

len_low = len(low_0_20_ts_dfs_A)
len_high = len(high_0_20_ts_dfs_A)

for i in range ((len_low//2)+1):
    try:
        j = i*2
        k = j+2
        # fname = 'LowNanp020comparisonA_TS_FS_normed'+str(j)
        # TS_and_FS(low_0_20_ts_dfs_A[j:k], low_0_20_fs_dfs_normed_A[j:k], low_0_20_labels_A[j:k], tak_freq=low_0_20_tak_A[j:k],norm=True, filename=fname)
        print('Cone Angle:'+str(low_0_20_cone_angle_A[j])+ ', pc4/3 ratio: '+str(low_0_20_ratioPc4_3_perp_A[j])+', Bomni:'+str(low_0_20_B_IMF_A[j])+' nT, Vsw:'+str(low_0_20_V_SW_A[j])+' km/s, Proton Density:'+str(low_0_20_n_P_A[j])+'/cc')
        print('Cone Angle:'+str(low_0_20_cone_angle_A[j+1])+ ', pc4/3 ratio: '+str(low_0_20_ratioPc4_3_perp_A[j+1])+', Bomni:'+str(low_0_20_B_IMF_A[j+1])+' nT, Vsw:'+str(low_0_20_V_SW_A[j+1])+' km/s, Proton Density:'+str(low_0_20_n_P_A[j+1])+'/cc')
    except:
        ''

for i in range ((len_high//2)+1):
    try:
        j = i*2
        k = j+2
        #fname = 'LowNanp020comparisonA_TS_FS_normed'+str(j)
        #TS_and_FS(low_0_20_ts_dfs_A[j:k], low_0_20_fs_dfs_normed_A[j:k], low_0_20_labels_A[j:k], tak_freq=low_0_20_tak_A[j:k],norm=True, filename=fname)
        print('Cone Angle:'+str(high_0_20_cone_angle_A[j])+ ', pc4/3 ratio: '+str(high_0_20_ratioPc4_3_perp_A[j])+', Bomni:'+str(high_0_20_B_IMF_A[j])+' nT, Vsw:'+str(high_0_20_V_SW_A[j])+' km/s, Proton Density:'+str(high_0_20_n_P_A[j])+'/cc')
        print('Cone Angle:'+str(high_0_20_cone_angle_A[j+1])+', pc4/3 ratio: '+str(high_0_20_ratioPc4_3_perp_A[j+1])+', Bomni:'+str(high_0_20_B_IMF_A[j+1])+' nT, Vsw:'+str(high_0_20_V_SW_A[j+1])+' km/s, Proton Density:'+str(high_0_20_n_P_A[j+1])+'/cc')
    except:
        ''
        
all_low_0_20_A = pd.concat(low_0_20_fs_dfs_normed_A)

count_hist_low_0_20_A, bin_edg = np.histogram(all_low_0_20_A['Frequency'], bins=2688)
total_t_hist_low_0_20_A, bin_edg = np.histogram(all_low_0_20_A['Frequency'], bins=2688, weights=all_low_0_20_A['Total Transverse Power'])
mean_trans_hist_low_0_20_A = total_t_hist_low_0_20_A/count_hist_low_0_20_A
total_c_hist_low_0_20_A, bin_edg = np.histogram(all_low_0_20_A['Frequency'], bins=2688, weights=all_low_0_20_A['Compressive Power'])
mean_comp_hist_low_0_20_A = total_c_hist_low_0_20_A/count_hist_low_0_20_A

all_high_0_20_A = pd.concat(high_0_20_fs_dfs_normed_A)

count_hist_high_0_20_A, bin_edg = np.histogram(all_high_0_20_A['Frequency'], bins=2688)
total_t_hist_high_0_20_A, bin_edg = np.histogram(all_high_0_20_A['Frequency'], bins=2688, weights=all_high_0_20_A['Total Transverse Power'])
mean_trans_hist_high_0_20_A = total_t_hist_high_0_20_A/count_hist_high_0_20_A
total_c_hist_high_0_20_A, bin_edg = np.histogram(all_high_0_20_A['Frequency'], bins=2688, weights=all_high_0_20_A['Compressive Power'])
mean_comp_hist_high_0_20_A = total_c_hist_high_0_20_A/count_hist_high_0_20_A

median_trans_hist_high_0_20_A, bin_edg,_ = scipy.stats.binned_statistic(all_high_0_20_A['Frequency'], all_high_0_20_A['Total Transverse Power'], statistic='median', bins=2688)
median_comp_hist_high_0_20_A, bin_edg,_ = scipy.stats.binned_statistic(all_high_0_20_A['Frequency'], all_high_0_20_A['Compressive Power'], statistic='median', bins=2688)

#bin midpoints:

midpoint_list = np.arange(2.083333335e-03, 11.2, 4.16666667e-03)

#find the mean proton and alpha predicted frequencies
mean_tak_0_20_A_low = np.mean(np.array(low_0_20_tak_A))
mean_tak_0_20_A_high = np.mean(np.array(high_0_20_tak_A))
median_tak_0_20_A_high = np.median(np.array(high_0_20_tak_A))

spectra_overlay_plot_comparison(low_0_20_fs_dfs_normed_A, high_0_20_fs_dfs_normed_A, CS_name='Case Study 1, 0-20°', trans_ave_1=mean_trans_hist_low_0_20_A, comp_ave_1=mean_comp_hist_low_0_20_A, trans_ave_2=median_trans_hist_high_0_20_A, comp_ave_2=median_comp_hist_high_0_20_A, midpoints=midpoint_list, tak_freq_1=mean_tak_0_20_A_low, tak_freq_2=median_tak_0_20_A_high)

print('Low Alpha Mean Pc4/Pc3 ratios - Transverse:', round(np.mean(low_0_20_ratioPc4_3_perp_A), 2), ', Compressive:', round(np.mean(low_0_20_ratioPc4_3_para_A), 2))
print('High Alpha Mean Pc4/Pc3 ratios - Transverse:', round(np.mean(high_0_20_ratioPc4_3_perp_A), 2), ', Compressive:', round(np.mean(high_0_20_ratioPc4_3_para_A), 2))

#low Nanp, 20-40deg. 
path = '/Users/roseatkinson/Documents/Extreme_Spectra/'

fs_df_list_low_40 = []
win_starts_low_40 = []
sc_no_low_40 = []
labels_low_40 = []
omni_means_low_40 = []
x_means_low_40 = []
y_means_low_40 = []
error_files_low_40 = []
takahashi_low_40 = [] 
mean_B_cluster_low_40 = []
mag_comp_low_40 = []
v_sw_means_low_40 = []
cone_angle_low_40 = []
proton_density_low_40 = []

for datetime, row in cl_fs_lowNaNp_20_40.iterrows():
    int_start = datetime.strftime("%Y-%m-%d %H:%M:%S")
    sc_no = row['sc']
    filename = 'FS_'+int_start+'_'+sc_no+'.csv'
    filepath = path+filename
    try:
        f_df = pd.read_csv(filepath,encoding='utf-8')
    except:
        error_files_low_40.append(filename)
    else:
        fs_df_list_low_40.append(f_df)
        win_starts_low_40.append(datetime)
        alpha_comp = row['SW Na/Np (mean)']
        alpha_comp = round(alpha_comp, 3)
        label = r'$N\alpha /Np$=' + str(alpha_comp)
        labels_low_40.append(label)
        sc_no_low_40.append(sc_no)
        omni_mean = row['IMF B (mean)']
        omni_means_low_40.append(round(omni_mean,1))
        x_gipm_mean = row['GIPM X (OMNI mean)']
        y_gipm_mean = row['GIPM Y (OMNI mean)']
        x_means_low_40.append(x_gipm_mean)
        y_means_low_40.append(y_gipm_mean)
        tak_freq = row['Takahashi Frequency, Hz']
        takahashi_low_40.append(tak_freq)
        bmean_bomni = row['B mean/Bomni']
        mag_comp_low_40.append(bmean_bomni)
        cluster_b_mean = bmean_bomni*omni_mean
        mean_B_cluster_low_40.append(cluster_b_mean)
        v_sw_means_low_40.append(round(row['SW V (mean)']))
        cone_angle_low_40.append(round(row['cone angle (mean)']))
        proton_density_low_40.append(round(row['SW Np (mean)'],1))           
        
ts_dfs_low_40 = []

for win_start, sc in zip(win_starts_low_40, sc_no_low_40):
    win_end = win_start + dt.timedelta(minutes=4)
    sc_dict = {'C1': cl_ts_1, 'C2': cl_ts_2, 'C3': cl_ts_3, 'C4': cl_ts_4}
    ts_df = sc_dict[sc]
    ts_df_filt = ts_df.loc[((ts_df.index >= win_start) & (ts_df.index < win_end))]
    ts_dfs_low_40.append(ts_df_filt)

fs_dfs_normed_low_40 = []
int_perp_Pc4_low_40 = []
int_perp_Pc3_low_40 = []
int_para_Pc4_low_40 = []
int_para_Pc3_low_40 = []

for df, B in zip(fs_df_list_low_40, mean_B_cluster_low_40):
    new_df = df.copy()
    new_df['Total Transverse Power'] = new_df['Total Transverse Power']/(B**2)
    new_df['Compressive Power'] = new_df['Compressive Power']/(B**2)
    fs_dfs_normed_low_40.append(new_df)
    P_perp_low, P_perp_high, P_para_low, P_para_high = split_band_power_int(new_df)
    int_perp_Pc4_low_40.append(P_perp_low)
    int_perp_Pc3_low_40.append(P_perp_high)
    int_para_Pc4_low_40.append(P_para_low)
    int_para_Pc3_low_40.append(P_para_high)

ratio_Pc4_Pc3_perp_low_40 = np.array(int_perp_Pc4_low_40)/np.array(int_perp_Pc3_low_40)
ratio_Pc4_Pc3_para_low_40 = np.array(int_para_Pc4_low_40)/np.array(int_para_Pc3_low_40)

#spectra_overlay_plot(fs_dfs_normed_low_40, 'low', '20-40°')

x_locs_low_20_40 = np.array(x_means_low_40)
y_locs_low_20_40 = np.array(y_means_low_40)

xy_pairs_low_20_40 = np.stack((x_locs_low_20_40,y_locs_low_20_40), axis=-1)

#location_mapping(xy_pairs_low_20_40, "20–40°", filename='low_20_40_cases')

#high Nanp, 20-40deg. 
path = '/Users/roseatkinson/Documents/Extreme_Spectra/'

fs_df_list_high_40 = []
win_starts_high_40 = []
sc_no_high_40 = []
labels_high_40 = []
omni_means_high_40 = []
x_means_high_40 = []
y_means_high_40 = []
error_files_high_40 = []
takahashi_high_40 = []
mean_B_cluster_high_40 = []
mag_comp_high_40 = []
v_sw_means_high_40 = []
cone_angle_high_40 = []
proton_density_high_40 = []

for datetime, row in cl_fs_highNaNp_20_40.iterrows():
    int_start = datetime.strftime("%Y-%m-%d %H:%M:%S")
    sc_no = row['sc']
    filename = 'FS_'+int_start+'_'+sc_no+'.csv'
    filepath = path+filename
    try:
        f_df = pd.read_csv(filepath,encoding='utf-8')
    except:
        error_files_high_40.append(filename)
    else:
        fs_df_list_high_40.append(f_df)
        win_starts_high_40.append(datetime)
        alpha_comp = row['SW Na/Np (mean)']
        alpha_comp = round(alpha_comp, 3)
        label = r'$N\alpha /Np$=' + str(alpha_comp)
        labels_high_40.append(label)
        sc_no_high_40.append(sc_no)
        omni_mean = row['IMF B (mean)']
        omni_means_high_40.append(round(omni_mean,1))
        x_gipm_mean = row['GIPM X (OMNI mean)']
        y_gipm_mean = row['GIPM Y (OMNI mean)']
        x_means_high_40.append(x_gipm_mean)
        y_means_high_40.append(y_gipm_mean)
        tak_freq = row['Takahashi Frequency, Hz']
        takahashi_high_40.append(tak_freq)
        bmean_bomni = row['B mean/Bomni']
        mag_comp_high_40.append(bmean_bomni)
        cluster_b_mean = bmean_bomni*omni_mean
        mean_B_cluster_high_40.append(cluster_b_mean)
        v_sw_means_high_40.append(round(row['SW V (mean)']))
        cone_angle_high_40.append(round(row['cone angle (mean)']))
        proton_density_high_40.append(round(row['SW Np (mean)'],1))
        
ts_dfs_high_40 = []

for win_start, sc in zip(win_starts_high_40, sc_no_high_40):
    win_end = win_start + dt.timedelta(minutes=4)
    sc_dict = {'C1': cl_ts_1, 'C2': cl_ts_2, 'C3': cl_ts_3, 'C4': cl_ts_4}
    ts_df = sc_dict[sc]
    ts_df_filt = ts_df.loc[((ts_df.index >= win_start) & (ts_df.index < win_end))]
    ts_dfs_high_40.append(ts_df_filt)

fs_dfs_normed_high_40 = []
int_perp_Pc4_high_40 = []
int_perp_Pc3_high_40 = []
int_para_Pc4_high_40 = []
int_para_Pc3_high_40 = []

for df, B in zip(fs_df_list_low_40, mean_B_cluster_low_40):
    new_df = df.copy()
    new_df['Total Transverse Power'] = new_df['Total Transverse Power']/(B**2)
    new_df['Compressive Power'] = new_df['Compressive Power']/(B**2)
    fs_dfs_normed_high_40.append(new_df)
    P_perp_low, P_perp_high, P_para_low, P_para_high = split_band_power_int(new_df)
    int_perp_Pc4_high_40.append(P_perp_low)
    int_perp_Pc3_high_40.append(P_perp_high)
    int_para_Pc4_high_40.append(P_para_low)
    int_para_Pc3_high_40.append(P_para_high)

ratio_Pc4_Pc3_perp_high_40 = np.array(int_perp_Pc4_high_40)/np.array(int_perp_Pc3_high_40)
ratio_Pc4_Pc3_para_high_40 = np.array(int_para_Pc4_high_40)/np.array(int_para_Pc3_high_40)

#spectra_overlay_plot(fs_dfs_normed_high_40, 'high', '20-40°')

x_locs_high_20_40 = np.array(x_means_high_40)
y_locs_high_20_40 = np.array(y_means_high_40)

xy_pairs_high_20_40 = np.stack((x_locs_high_20_40,y_locs_high_20_40), axis=-1)

#location_mapping(xy_pairs_high_20_40, "20–40°", filename='high_20_40_cases')
#location_mapping(xy_pairs_low_20_40, "20–40°",  alt_locs=xy_pairs_high_20_40, filename='low_high_20_40_cases')

#locations:
#compare 14<X<17, -2<Y<0

#new set of lists

low_20_40_locations_A = []
low_20_40_ts_dfs_A = []
low_20_40_fs_dfs_normed_A = []
low_20_40_labels_A = []
low_20_40_tak_A = []
low_20_40_ratioPc4_3_perp_A = []
low_20_40_ratioPc4_3_para_A = []
low_20_40_V_SW_A = []
low_20_40_cone_angle_A = []
low_20_40_n_P_A = []
low_20_40_B_IMF_A = []

for ind in [20,21,22,23,29,39,48]:
    low_20_40_locations_A.append(xy_pairs_low_20_40[ind])
    low_20_40_ts_dfs_A.append(ts_dfs_low_40[ind])
    low_20_40_labels_A.append(labels_low_40[ind])
    low_20_40_tak_A.append(takahashi_low_40[ind])
    low_20_40_fs_dfs_normed_A.append(fs_dfs_normed_low_40[ind])
    low_20_40_ratioPc4_3_perp_A.append(ratio_Pc4_Pc3_perp_low_40[ind])
    low_20_40_ratioPc4_3_para_A.append(ratio_Pc4_Pc3_para_low_40[ind])
    low_20_40_V_SW_A.append(v_sw_means_low_40[ind])
    low_20_40_cone_angle_A.append(cone_angle_low_40[ind])
    low_20_40_n_P_A.append(proton_density_low_40[ind])
    low_20_40_B_IMF_A.append(omni_means_low_40[ind])

high_20_40_locations_A = []
high_20_40_ts_dfs_A = []
high_20_40_fs_dfs_normed_A = []
high_20_40_labels_A = []
high_20_40_tak_A = []
high_20_40_ratioPc4_3_perp_A = []
high_20_40_ratioPc4_3_para_A = []
high_20_40_V_SW_A = []
high_20_40_cone_angle_A = []
high_20_40_n_P_A = []
high_20_40_B_IMF_A = []

for ind in [9,16,20,21,22,27,28,29]:
    high_20_40_locations_A.append(xy_pairs_high_20_40[ind])
    high_20_40_ts_dfs_A.append(ts_dfs_high_40[ind])
    high_20_40_labels_A.append(labels_high_40[ind])
    high_20_40_tak_A.append(takahashi_high_40[ind])
    high_20_40_fs_dfs_normed_A.append(fs_dfs_normed_high_40[ind])
    high_20_40_ratioPc4_3_perp_A.append(ratio_Pc4_Pc3_perp_high_40[ind])
    high_20_40_ratioPc4_3_para_A.append(ratio_Pc4_Pc3_para_high_40[ind])
    high_20_40_V_SW_A.append(v_sw_means_high_40[ind])
    high_20_40_cone_angle_A.append(cone_angle_high_40[ind])
    high_20_40_n_P_A.append(proton_density_high_40[ind])
    high_20_40_B_IMF_A.append(omni_means_high_40[ind])

location_mapping(low_20_40_locations_A, "20–40°",  alt_locs=high_20_40_locations_A, filename='low_high_20_40_cases_A')

len_low = len(low_20_40_ts_dfs_A)
len_high = len(high_20_40_ts_dfs_A)

for i in range ((len_low//2)+1):
    try:
        j = i*2
        k = j+2
        #fname = 'LowNanp2040comparisonA_TS_FS_normed'+str(j)
        #TS_and_FS(low_20_40_ts_dfs_A[j:k], low_20_40_fs_dfs_normed_A[j:k], low_20_40_labels_A[j:k], tak_freq=low_20_40_tak_A[j:k],norm=True, filename=fname)
        print('Cone Angle:'+str(low_20_40_cone_angle_A[j])+', pc4/3 ratio: '+str(low_20_40_ratioPc4_3_perp_A[j])+', Bomni:'+str(low_20_40_B_IMF_A[j])+' nT, Vsw:'+str(low_20_40_V_SW_A[j])+' km/s, Proton Density:'+str(high_20_40_n_P_A[j])+'/cc')
        print('Cone Angle:'+str(low_20_40_cone_angle_A[j+1])+', pc4/3 ratio: '+str(low_20_40_ratioPc4_3_perp_A[j+1])+', Bomni:'+str(low_20_40_B_IMF_A[j+1])+' nT, Vsw:'+str(low_20_40_V_SW_A[j+1])+' km/s, Proton Density:'+str(high_20_40_n_P_A[j+1])+'/cc')
    except:
        ''
        
all_low_20_40_A = pd.concat(low_20_40_fs_dfs_normed_A)

count_hist_low_20_40_A, bin_edg = np.histogram(all_low_20_40_A['Frequency'], bins=2688)
total_t_hist_low_20_40_A, bin_edg = np.histogram(all_low_20_40_A['Frequency'], bins=2688, weights=all_low_20_40_A['Total Transverse Power'])
mean_trans_hist_low_20_40_A = total_t_hist_low_20_40_A/count_hist_low_20_40_A
total_c_hist_low_20_40_A, bin_edg = np.histogram(all_low_20_40_A['Frequency'], bins=2688, weights=all_low_20_40_A['Compressive Power'])
mean_comp_hist_low_20_40_A = total_c_hist_low_20_40_A/count_hist_low_20_40_A

median_trans_hist_low_20_40_A, bin_edg,_ = scipy.stats.binned_statistic(all_low_20_40_A['Frequency'], all_low_20_40_A['Total Transverse Power'], statistic='median', bins=2688)
median_comp_hist_low_20_40_A, bin_edg,_ = scipy.stats.binned_statistic(all_low_20_40_A['Frequency'], all_low_20_40_A['Compressive Power'], statistic='median', bins=2688)

all_high_20_40_A = pd.concat(high_20_40_fs_dfs_normed_A)

count_hist_high_20_40_A, bin_edg = np.histogram(all_high_20_40_A['Frequency'], bins=2688)
total_t_hist_high_20_40_A, bin_edg = np.histogram(all_high_20_40_A['Frequency'], bins=2688, weights=all_high_20_40_A['Total Transverse Power'])
mean_trans_hist_high_20_40_A = total_t_hist_high_20_40_A/count_hist_high_20_40_A
total_c_hist_high_20_40_A, bin_edg = np.histogram(all_high_20_40_A['Frequency'], bins=2688, weights=all_high_20_40_A['Compressive Power'])
mean_comp_hist_high_20_40_A = total_c_hist_high_20_40_A/count_hist_high_20_40_A

median_trans_hist_high_20_40_A, bin_edg,_ = scipy.stats.binned_statistic(all_high_20_40_A['Frequency'], all_high_20_40_A['Total Transverse Power'], statistic='median', bins=2688)
median_comp_hist_high_20_40_A, bin_edg,_ = scipy.stats.binned_statistic(all_high_20_40_A['Frequency'], all_high_20_40_A['Compressive Power'], statistic='median', bins=2688)

#bin midpoints:

midpoint_list = np.arange(2.083333335e-03, 11.2, 4.16666667e-03)

#spectra_overlay_plot(low_20_40_fs_dfs_normed_A, 'low', '20-40°', trans_ave=mean_trans_hist_low_20_40_a, comp_ave=mean_comp_hist_low_20_40_a, midpoints=midpoint_list)
#spectra_overlay_plot(high_20_40_fs_dfs_normed_A, 'high', '20-40°', trans_ave=mean_trans_hist_high_20_40_a, comp_ave=mean_comp_hist_high_20_40_a, midpoints=midpoint_list)

#find the mean proton and alpha predicted frequencies
mean_tak_20_40_A_low = np.median(np.array(low_20_40_tak_A))
mean_tak_20_40_A_high = np.median(np.array(high_20_40_tak_A))
median_tak_20_40_A_low = np.median(np.array(low_20_40_tak_A))
median_tak_20_40_A_high = np.median(np.array(high_20_40_tak_A))

spectra_overlay_plot_comparison(low_20_40_fs_dfs_normed_A, high_20_40_fs_dfs_normed_A, CS_name='Multi Case Study 2A, 20-40°', trans_ave_1=median_trans_hist_low_20_40_A, comp_ave_1=median_comp_hist_low_20_40_A, trans_ave_2=median_trans_hist_high_20_40_A, comp_ave_2=median_comp_hist_high_20_40_A, midpoints=midpoint_list, tak_freq_1=median_tak_20_40_A_low, tak_freq_2=median_tak_20_40_A_high)

print('Low Alpha Mean Pc4/Pc3 ratios - Transverse:', round(np.mean(low_20_40_ratioPc4_3_perp_A), 2), ', Compressive:', round(np.mean(low_20_40_ratioPc4_3_para_A), 2))
print('High Alpha Mean Pc4/Pc3 ratios - Transverse:', round(np.mean(high_20_40_ratioPc4_3_perp_A), 2), ', Compressive:', round(np.mean(high_20_40_ratioPc4_3_para_A), 2))

for i in range ((len_high//2)+1):
    try:
        j = i*2
        k = j+2
        # fname = 'highNanp2040comparisonA_TS_FS_normed'+str(j)
        # TS_and_FS(high_20_40_ts_dfs_A[j:k], high_20_40_fs_dfs_normed_A[j:k], high_20_40_labels_A[j:k], tak_freq=high_20_40_tak_A[j:k],norm=True, filename=fname)
        print('Cone Angle:'+str(high_20_40_cone_angle_A[j])+', pc4/3 ratio: '+str(high_20_40_ratioPc4_3_perp_A[j])+ ', Bomni:'+str(high_20_40_B_IMF_A[j])+' nT, Vsw:'+str(high_20_40_V_SW_A[j])+' km/s, Proton Density:'+str(high_20_40_n_P_A[j])+'/cc')
        print('Cone Angle:'+str(high_20_40_cone_angle_A[j+1])+', pc4/3 ratio: '+str(high_20_40_ratioPc4_3_perp_A[j+1])+ ', Bomni:'+str(high_20_40_B_IMF_A[j+1])+' nT, Vsw:'+str(high_20_40_V_SW_A[j+1])+' km/s, Proton Density:'+str(high_20_40_n_P_A[j+1])+'/cc')
    except:
        print(j,k)

# B: 14<X<17, -4<Y<-2, 
#print(xy_pairs_low_20_40[19], xy_pairs_low_20_40[24], xy_pairs_low_20_40[27], xy_pairs_low_20_40[30:34], xy_pairs_low_20_40[49])

#new set of lists

low_20_40_locations_B = []
low_20_40_ts_dfs_B = []
low_20_40_fs_dfs_normed_B = []
low_20_40_labels_B = []
low_20_40_tak_B = []
low_20_40_ratioPc4_3_perp_B = []
low_20_40_ratioPc4_3_para_B = []
low_20_40_V_SW_B = []
low_20_40_cone_angle_B = []
low_20_40_n_P_B = []
low_20_40_B_IMF_B = []

for ind in [19,24,27,30,31,32,33,49]:
    low_20_40_locations_B.append(xy_pairs_low_20_40[ind])
    low_20_40_ts_dfs_B.append(ts_dfs_low_40[ind])
    low_20_40_labels_B.append(labels_low_40[ind])
    low_20_40_tak_B.append(takahashi_low_40[ind])
    low_20_40_fs_dfs_normed_B.append(fs_dfs_normed_low_40[ind])
    low_20_40_ratioPc4_3_perp_B.append(ratio_Pc4_Pc3_perp_low_40[ind])
    low_20_40_ratioPc4_3_para_B.append(ratio_Pc4_Pc3_para_low_40[ind])
    low_20_40_V_SW_B.append(v_sw_means_low_40[ind])
    low_20_40_cone_angle_B.append(cone_angle_low_40[ind])
    low_20_40_n_P_B.append(proton_density_low_40[ind])
    low_20_40_B_IMF_B.append(omni_means_low_40[ind])

high_20_40_locations_B = []
high_20_40_ts_dfs_B = []
high_20_40_fs_dfs_normed_B = []
high_20_40_labels_B = []
high_20_40_tak_B = []
high_20_40_ratioPc4_3_perp_B = []
high_20_40_ratioPc4_3_para_B = []
high_20_40_V_SW_B = []
high_20_40_cone_angle_B = []
high_20_40_n_P_B = []
high_20_40_B_IMF_B = []

for ind in [10,11,12,14,15,23,24,25,26,30,31,32,33]:
    high_20_40_locations_B.append(xy_pairs_high_20_40[ind])
    high_20_40_ts_dfs_B.append(ts_dfs_high_40[ind])
    high_20_40_labels_B.append(labels_high_40[ind])
    high_20_40_tak_B.append(takahashi_high_40[ind])
    high_20_40_fs_dfs_normed_B.append(fs_dfs_normed_high_40[ind])
    high_20_40_ratioPc4_3_perp_B.append(ratio_Pc4_Pc3_perp_high_40[ind])
    high_20_40_ratioPc4_3_para_B.append(ratio_Pc4_Pc3_para_high_40[ind])
    high_20_40_V_SW_B.append(v_sw_means_high_40[ind])
    high_20_40_cone_angle_B.append(cone_angle_high_40[ind])
    high_20_40_n_P_B.append(proton_density_high_40[ind])
    high_20_40_B_IMF_B.append(omni_means_high_40[ind])
    
location_mapping(low_20_40_locations_B, "20–40°",  alt_locs=high_20_40_locations_B, filename='low_high_20_40_cases_B')

len_low = len(low_20_40_ts_dfs_B)
len_high = len(high_20_40_ts_dfs_B)

for i in range ((len_low//2)+1):
    try:
        j = i*2
        k = j+2
        #fname = 'LowNanp2040comparisonB_TS_FS_normed'+str(j)
        #TS_and_FS(low_20_40_ts_dfs_B[j:k], low_20_40_fs_dfs_normed_B[j:k], low_20_40_labels_B[j:k], tak_freq=low_20_40_tak_B[j:k],norm=True, filename=fname)
    except:
        'run out!'

all_low_20_40_B = pd.concat(low_20_40_fs_dfs_normed_B)

count_hist_low_20_40_B, bin_edg = np.histogram(all_low_20_40_B['Frequency'], bins=2688)
total_t_hist_low_20_40_B, bin_edg = np.histogram(all_low_20_40_B['Frequency'], bins=2688, weights=all_low_20_40_B['Total Transverse Power'])
mean_trans_hist_low_20_40_B = total_t_hist_low_20_40_B/count_hist_low_20_40_B
total_c_hist_low_20_40_B, bin_edg = np.histogram(all_low_20_40_B['Frequency'], bins=2688, weights=all_low_20_40_B['Compressive Power'])
mean_comp_hist_low_20_40_B = total_c_hist_low_20_40_B/count_hist_low_20_40_B

median_trans_hist_low_20_40_B, bin_edg,_ = scipy.stats.binned_statistic(all_low_20_40_B['Frequency'], all_low_20_40_B['Total Transverse Power'], statistic='median', bins=2688)
median_comp_hist_low_20_40_B, bin_edg,_ = scipy.stats.binned_statistic(all_low_20_40_B['Frequency'], all_low_20_40_B['Compressive Power'], statistic='median', bins=2688)

all_high_20_40_B = pd.concat(high_20_40_fs_dfs_normed_B)

count_hist_high_20_40_B, bin_edg = np.histogram(all_high_20_40_B['Frequency'], bins=2688)
total_t_hist_high_20_40_B, bin_edg = np.histogram(all_high_20_40_B['Frequency'], bins=2688, weights=all_high_20_40_B['Total Transverse Power'])
mean_trans_hist_high_20_40_B = total_t_hist_high_20_40_B/count_hist_high_20_40_B
total_c_hist_high_20_40_B, bin_edg = np.histogram(all_high_20_40_B['Frequency'], bins=2688, weights=all_high_20_40_B['Compressive Power'])
mean_comp_hist_high_20_40_B = total_c_hist_high_20_40_B/count_hist_high_20_40_B

median_trans_hist_high_20_40_B, bin_edg,_ = scipy.stats.binned_statistic(all_high_20_40_B['Frequency'], all_high_20_40_B['Total Transverse Power'], statistic='median', bins=2688)
median_comp_hist_high_20_40_B, bin_edg,_ = scipy.stats.binned_statistic(all_high_20_40_B['Frequency'], all_high_20_40_B['Compressive Power'], statistic='median', bins=2688)

#bin midpoints:

midpoint_list = np.arange(2.083333335e-03, 11.2, 4.16666667e-03)

#spectra_overlay_plot(low_20_40_fs_dfs_normed_B, 'low', '20-40°', trans_ave=mean_trans_hist_low_20_40_a, comp_ave=mean_comp_hist_low_20_40_a, midpoints=midpoint_list)
#spectra_overlay_plot(high_20_40_fs_dfs_normed_B, 'high', '20-40°', trans_ave=mean_trans_hist_high_20_40_a, comp_ave=mean_comp_hist_high_20_40_a, midpoints=midpoint_list)

#find the mean proton and alpha predicted frequencies
mean_tak_20_40_B_low = np.median(np.array(low_20_40_tak_B))
mean_tak_20_40_B_high = np.median(np.array(high_20_40_tak_B))
median_tak_20_40_B_low = np.median(np.array(low_20_40_tak_B))
median_tak_20_40_B_high = np.median(np.array(high_20_40_tak_B))

spectra_overlay_plot_comparison(low_20_40_fs_dfs_normed_B, high_20_40_fs_dfs_normed_B, CS_name='Case Study 2B, 20-40°', trans_ave_1=median_trans_hist_low_20_40_B, comp_ave_1=median_comp_hist_low_20_40_B, trans_ave_2=median_trans_hist_high_20_40_B, comp_ave_2=median_comp_hist_high_20_40_B, midpoints=midpoint_list, tak_freq_1=median_tak_20_40_B_low, tak_freq_2=median_tak_20_40_B_high)


len_low = len(low_20_40_ts_dfs_B)
len_high = len(high_20_40_ts_dfs_B)

for i in range ((len_low//2)+1):
    try:
        j = i*2
        k = j+2
        print('Cone Angle:'+str(low_20_40_cone_angle_B[j])+', pc4/3 ratio: '+str(low_20_40_ratioPc4_3_perp_B[j])+', Bomni:'+str(low_20_40_B_IMF_B[j])+' nT, Vsw:'+str(low_20_40_V_SW_B[j])+' km/s, Proton Density:'+str(low_20_40_n_P_B[j])+'/cc')
        print('Cone Angle:'+str(low_20_40_cone_angle_B[j+1])+', pc4/3 ratio: '+str(low_20_40_ratioPc4_3_perp_B[j+1])+', Bomni:'+str(low_20_40_B_IMF_B[j+1])+' nT, Vsw:'+str(low_20_40_V_SW_B[j+1])+' km/s, Proton Density:'+str(low_20_40_n_P_B[j+1])+'/cc')
    except:
        ''
        
#spectra_overlay_plot_comparison(low_20_40_fs_dfs_normed_B, high_20_40_fs_dfs_normed_B, '20-40°', fname='locB', trans_ave_1=mean_trans_hist_low_20_40_B, comp_ave_1=mean_comp_hist_low_20_40_B, trans_ave_2=mean_trans_hist_high_20_40_B, comp_ave_2=mean_comp_hist_high_20_40_B, midpoints=midpoint_list, tak_freq_1=mean_tak_20_40_B_low, tak_freq_2=mean_tak_20_40_B_high)

print('Low Alpha Mean Pc4/Pc3 ratios - Transverse:', round(np.mean(low_20_40_ratioPc4_3_perp_B),2), ', Compressive:',round(np.mean(low_20_40_ratioPc4_3_para_B),2))
print('High Alpha Mean Pc4/Pc3 ratios - Transverse:', round(np.mean(high_20_40_ratioPc4_3_perp_B),2), ', Compressive:', round(np.mean(high_20_40_ratioPc4_3_para_B),2))
