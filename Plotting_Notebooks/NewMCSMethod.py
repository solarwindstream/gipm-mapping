import pandas as pd
import numpy as np
import datetime as dt
import glob

import matplotlib 
from matplotlib import pyplot as plt
from spectra_modules import spectra_plot, TS_plot, spectra_overlay_plot, TS_and_FS, split_band_power_int, spectra_overlay_plot_comparison
from gipm_plotting_modules import location_mapping
from heatmap_modules import compute_hists2d
import scipy

##load integrated Cluster CSVs

path = "/Users/roseatkinson/Documents/Cluster_Int_CSVs/**"

cl_file_list = {'C1': [], 'C2': [], 'C3': [], 'C4': []}

for path in glob.glob(path, recursive=True):
    if '.csv' in path:
        if 'C1' in path:
            cl_file_list['C1'].append(path)
        if 'C2' in path:
            cl_file_list['C2'].append(path)
        if 'C3' in path:
            cl_file_list['C3'].append(path)
        if 'C4' in path:
            cl_file_list['C4'].append(path)
    
cl_dfs = []
scs = ['C1', 'C2', 'C3', 'C4']

for sc in scs:
    for file in cl_file_list[sc]:
        df = pd.read_csv(file,encoding='utf-8')
        df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
        df.set_index('datetime', inplace = True)
        df['sc']=sc
        cl_dfs.append(df)

cl_power_all = pd.concat(cl_dfs)

print("dfs loaded in")

#Filtering: remove times when OMNI source sc was too far from X-line, & current sheets, & filter down to -5<Zgipm<5 "plane"
cl_filtered = cl_power_all.loc[(cl_power_all['OMNI Dist from X line (mean)'] < 70) & (cl_power_all['Max IMF Deviation'] < 60) & (cl_power_all['GIPM Z (OMNI mean)'] < 5) & (cl_power_all['GIPM Z (OMNI mean)'] > -5)]
#filtered for 3<B<6
cl_B_filtered = cl_filtered.loc[(cl_filtered['IMF B (mean)'] > 3) & (cl_filtered['IMF B (mean)'] < 6)]
#filter for likely ULF
cl_B_filtered_ULF = cl_B_filtered.loc[(cl_B_filtered['ULF Band Normalised Transverse Power'] > 0.05)]
# & filter for cone angle class
cl_filtered_0_20 = cl_B_filtered_ULF.loc[(cl_B_filtered_ULF['cone angle (mean)'] >= 0) & (cl_B_filtered_ULF['cone angle (mean)'] < 20)]
cl_filtered_20_40 = cl_B_filtered_ULF.loc[(cl_B_filtered_ULF['cone angle (mean)'] >= 20) & (cl_B_filtered_ULF['cone angle (mean)'] < 40)]

##separately, load time series CSVs for extreme alpha cases
ex_cl_df_list = {'C1': [], 'C2': [], 'C3': [], 'C4': []}

path = "/Users/roseatkinson/Documents/Cluster_Time_Series/**"

for path in glob.glob(path, recursive=True):
    if '.csv' in path:
        df = pd.read_csv(path,encoding='utf-8')
        df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
        df.set_index('datetime', inplace = True)
        if 'C1' in path:
            ex_cl_df_list['C1'].append(df)
        if 'C2' in path:
            ex_cl_df_list['C2'].append(df)
        if 'C3' in path:
            ex_cl_df_list['C3'].append(df)
        if 'C4' in path:
            ex_cl_df_list['C4'].append(df)

#produce one time series df per spacecraft sorted by time
cl_ts = {}

for sc in scs:
    cl_ts[sc]=pd.concat(ex_cl_df_list[sc])
    cl_ts[sc].sort_index(inplace=True)

#print(cl_ts['C1'].index)

#for finding likely foreshock areas, no longer needed but for documentation purposes
#from gipm_plotting_modules import single_map
# heatmaps_dict_new = {'0-20°':{},'20-40°':{}}

# heatmaps_dict_new['0-20°'] = compute_hists2d(cl_filtered_0_20, 'mean', 'Normalised Transverse Power', obs_min=50)
# heatmaps_dict_new['20-40°'] = compute_hists2d(cl_filtered_20_40, 'mean', 'Normalised Transverse Power', obs_min=50)

#note that ion and tangent lines are for 30-45 case
#single_map(heatmaps_dict_new['0-20°']['Normalised Transverse Power'], "0–20°",'Power')
#single_map(heatmaps_dict_new['20-40°']['Normalised Transverse Power'], "20–40°",'Power')

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

fs_df_list_0_20 = []
fs_df_list_20_40 = []

for x in range(9,16):
    x_end = x + 1
    y_bounds = y_lims_0_20[str(x)]
    df = cl_filtered_0_20.loc[(cl_filtered_0_20['GIPM X (OMNI mean)'] >= x) & (cl_filtered_0_20['GIPM X (OMNI mean)'] < x_end) & (cl_filtered_0_20['GIPM Y (OMNI mean)'] >= y_bounds[0]) & (cl_filtered_0_20['GIPM Y (OMNI mean)'] < y_bounds[1])]
    fs_df_list_0_20.append(df)

for x in range(5,17):
    x_end = x + 1
    y_bounds = y_lims_20_40[str(x)]
    df = cl_filtered_20_40_ULF.loc[(cl_filtered_20_40_ULF['GIPM X (OMNI mean)'] >= x) & (cl_filtered_20_40_ULF['GIPM X (OMNI mean)'] < x_end) & (cl_filtered_20_40_ULF['GIPM Y (OMNI mean)'] >= y_bounds[0]) & (cl_filtered_20_40_ULF['GIPM Y (OMNI mean)'] < y_bounds[1])]
    fs_df_list_20_40.append(df)

#Background start times and mean background conditions for all possible foreshock
fs_dfs_0_20 = pd.concat(fs_df_list_0_20)
fs_dfs_20_40 = pd.concat(fs_df_list_20_40)

high_alpha=0.08
low_alpha=0.005

cl_intervals_foreshock = {'0_20':{},'20_40': {}}

cl_intervals_foreshock['0_20']['High'] = fs_dfs_0_20.loc[(fs_dfs_0_20['SW Na/Np (mean)'] > high_alpha)]
cl_intervals_foreshock['0_20']['Low'] = fs_dfs_0_20.loc[(fs_dfs_0_20['SW Na/Np (mean)'] < low_alpha)]
cl_intervals_foreshock['20_40']['High'] = fs_dfs_20_40.loc[(fs_dfs_20_40['SW Na/Np (mean)'] > high_alpha)]
cl_intervals_foreshock['20_40']['Low'] = fs_dfs_20_40.loc[(fs_dfs_20_40['SW Na/Np (mean)'] < low_alpha)]

#Now for each interval pair up the Fourier spectrum and corresponding slice of the Cluster time series
#as well as the mean background conditions & location for that interval

path = '/Users/roseatkinson/Documents/Extreme_Spectra/'

#Dictionary for NaNp/cone angle bins

CA_NANP_binning = {}

errorfiles = []


#pairings of group name and interval list dfs to consult

for ca, alpha_dict in cl_intervals_foreshock.items():
    CA_NANP_binning[ca]={}
        for alpha, df in alpha_dict.items():
            CA_NANP_binning[ca][alpha]={}
            #now 
            
    fs_df_list, fs_normed_list, win_starts, sc_no_list, ts_df_list, nanp, labels, omni_means, x_means, y_means, takahashi, mean_B_cluster, mag_comp, v_sw_means, cone_angle, proton_density, tak_trans_diff, tak_comp_diff, peak_t_freq, peak_c_freq, ratio_Pc4_Pc3_trans, ratio_Pc4_Pc3_comp  =([] for i in range(22))
    
    for datetime, row in df.iterrows():
        int_start = datetime.strftime("%Y-%m-%d %H:%M:%S")
        sc_no = row['sc']
        filename = 'FS_'+int_start+'_'+sc_no+'.csv'
        filepath = path+filename
        try:
            f_df = pd.read_csv(filepath,encoding='utf-8')
        except:
            errorfiles.append(filename)
        else:
            fs_df_list.append(f_df)
            win_starts.append(datetime)
            #select the four-minute window of raw time series data
            win_end = datetime + dt.timedelta(minutes=4)
            ts_df = sc_dict[sc_no]
            ts_df_filt = ts_df.loc[((ts_df.index >= datetime) & (ts_df.index < win_end))]
            ts_df_list.append(ts_df_filt)
            #find upstream & local parameter values for interval 
            alpha_comp = row['SW Na/Np (mean)']
            nanp.append(alpha_comp)
            label = r'$N\alpha /Np$=' + str(round(alpha_comp, 3))
            labels.append(label)
            sc_no_list.append(sc_no)
            omni_mean = row['IMF B (mean)']
            omni_means.append(round(omni_mean,1))
            x_gipm_mean = row['GIPM X (OMNI mean)']
            y_gipm_mean = row['GIPM Y (OMNI mean)']
            tak_freq = row['Takahashi Frequency, Hz']
            takahashi.append(tak_freq)
            x_means.append(x_gipm_mean)
            y_means.append(y_gipm_mean)
            bmean_bomni = row['B mean/Bomni']
            mag_comp.append(bmean_bomni)
            cluster_b_mean = bmean_bomni*omni_mean
            mean_B_cluster.append(cluster_b_mean)
            v_sw_means.append(round(row['SW V (mean)']))
            cone_angle.append(round(row['cone angle (mean)']))
            proton_density.append(round(row['SW Np (mean)'], 1)) 
            #transform units to mHz from Hz, & reverse diff from Ftak-Fpeak to Fpeak-Ftak:
            tak_trans_diff.append(-1000*row['Takahashi Transverse Difference'])
            tak_comp_diff.append(-1000*row['Takahashi Compressive Difference'])
            peak_t_freq.append(1000*row['Peak Transverse Frequency'])
            peak_c_freq.append(1000*row['Peak Compressive Frequency'])
            #Also generate normalised power spectrum!
            new_df = f_df.copy()
            new_df['Total Transverse Power'] = new_df['Total Transverse Power']/(cluster_b_mean**2)
            new_df['Compressive Power'] = new_df['Compressive Power']/(cluster_b_mean**2)
            fs_normed_list.append(new_df)
            P_perp_low, P_perp_high, P_para_low, P_para_high = split_band_power_int(new_df)
            ratio_trans = P_perp_low/P_perp_high
            ratio_comp = P_para_low/P_para_high
            ratio_Pc4_Pc3_trans.append(ratio_trans)
            ratio_Pc4_Pc3_comp.append(ratio_Pc4_Pc3_comp)
            
    x_locs = np.array(x_means)
    y_locs = np.array(y_means)

    xy_pairs = np.stack((x_locs,y_locs), axis=-1)

    CA_NANP_binning[ref]['TS List']=ts_df_list
    CA_NANP_binning[ref]['FS List']=fs_df_list
    CA_NANP_binning[ref]['Normed FS List']=fs_normed_list
    CA_NANP_binning[ref]['Win Starts']=win_starts
    CA_NANP_binning[ref]['Sc No.']= sc_no_list
    CA_NANP_binning[ref]['NaNp']=nanp
    CA_NANP_binning[ref]['NaNp Labels']=labels
    CA_NANP_binning[ref]['OMNI means']=omni_means
    CA_NANP_binning[ref]['GIPM X mean']=x_means
    CA_NANP_binning[ref]['GIPM Y mean']=y_means
    CA_NANP_binning[ref]['Takahashi Frequency']=takahashi
    CA_NANP_binning[ref]['Cluster Mean B']=mean_B_cluster
    CA_NANP_binning[ref]['Bcluster/Bomni']=mag_comp
    CA_NANP_binning[ref]['Mean V_SW']=v_sw_means
    CA_NANP_binning[ref]['Cone Angle']=cone_angle
    CA_NANP_binning[ref]['Np']=proton_density
    CA_NANP_binning[ref]['Takahashi Transverse Difference']=tak_trans_diff
    CA_NANP_binning[ref]['Takahashi Compressive Difference']=tak_comp_diff
    CA_NANP_binning[ref]['Peak Transverse Frequency']=peak_t_freq
    CA_NANP_binning[ref]['Peak Compressive Frequency']=peak_c_freq
    CA_NANP_binning[ref]['Transverse Pc4/Pc3 Ratio']=ratio_Pc4_Pc3_trans
    CA_NANP_binning[ref]['Compressive Pc4/Pc3 Ratio']=ratio_Pc4_Pc3_comp
    CA_NANP_binning[ref]['XY mean locations']=xy_pairs

all_low_mean_b = []
all_high_mean_b = []
all_low_tak = []
all_high_tak = []

for label_dict in MCS_Dicts_Labels:
    all_low_mean_b.extend(MCS_Dicts[label_dict]['Low']['OMNI means'])
    all_high_mean_b.extend(MCS_Dicts[label_dict]['High']['OMNI means'])
    all_low_tak.extend(MCS_Dicts[label_dict]['Low']['Takahashi Frequency'])
    all_high_tak.extend(MCS_Dicts[label_dict]['High']['Takahashi Frequency'])
    
    low_a_mean_B = np.mean(np.array(MCS_Dicts[label_dict]['Low']['OMNI means']))
    high_a_mean_B = np.mean(np.array(MCS_Dicts[label_dict]['High']['OMNI means']))
    low_a_tak_mean = np.mean(np.array(MCS_Dicts[label_dict]['Low']['Takahashi Frequency']))
    high_a_tak_mean = np.mean(np.array(MCS_Dicts[label_dict]['High']['Takahashi Frequency']))

    #print(f"{label_dict}: High-alpha B, Low-alpha B: {high_a_mean_B} {low_a_mean_B}. High-alpha Ftak, Low-alpha Ftak: {high_a_tak_mean}{low_a_tak_mean}")

overall_low_b = np.mean(np.array(all_low_mean_b))
overall_high_b = np.mean(np.array(all_high_mean_b))
overall_low_tak = np.mean(np.array(all_low_tak))
overall_high_tak = np.mean(np.array(all_high_tak))

print(f"Overall: High-alpha B, Low-alpha B: {overall_high_b}  {overall_low_b}.")
      
print(f"Overall: High-alpha Ftak, Low-alpha Ftak: {overall_high_tak*1000}  {overall_low_tak*1000}")

#find the mean & median of the spectra

bin_list = np.arange(0.0041627950195131/2, 11.2, 0.0041627950195131)
midpoint_list= np.arange(0.0041627950195131, 11.2, 0.0041627950195131)

MCS_FS_Dicts = {}

#mean

for label in MCS_Dicts_Labels:
    
    MCS_FS_Dicts[label]={}
    MCS_FS_Dicts[label]['Low']={}
    MCS_FS_Dicts[label]['High']={}

    #Mean
    
    all_fs_low = pd.concat(MCS_Dicts[label]['Low']['Normed FS List'])
    count_hist_low, bin_edg = np.histogram(all_fs_low['Frequency'], bins=bin_list)
    total_t_hist_low, bin_edg = np.histogram(all_fs_low['Frequency'], bins=bin_list, weights=all_fs_low['Total Transverse Power'])
    MCS_FS_Dicts[label]['Low']['Mean Transverse Hist'] = total_t_hist_low/count_hist_low
    total_c_hist_low, bin_edg = np.histogram(all_fs_low['Frequency'], bins=bin_list, weights=all_fs_low['Compressive Power'])
    MCS_FS_Dicts[label]['Low']['Mean Compressive Hist'] = total_c_hist_low/count_hist_low

    all_fs_high = pd.concat(MCS_Dicts[label]['High']['Normed FS List'])
    count_hist_high, bin_edg = np.histogram(all_fs_high['Frequency'], bins=bin_list)
    total_t_hist_high, bin_edg = np.histogram(all_fs_high['Frequency'], bins=bin_list, weights=all_fs_high['Total Transverse Power'])
    MCS_FS_Dicts[label]['High']['Mean Transverse Hist'] = total_t_hist_high/count_hist_high
    total_c_hist_high, bin_edg = np.histogram(all_fs_high['Frequency'], bins=bin_list, weights=all_fs_high['Compressive Power'])
    MCS_FS_Dicts[label]['High']['Mean Compressive Hist'] = total_c_hist_high/count_hist_high

    #Median
    MCS_FS_Dicts[label]['Low']['Median Transverse Hist'], bin_edg,_ = scipy.stats.binned_statistic(all_fs_low['Frequency'], all_fs_low['Total Transverse Power'], statistic='median', bins=bin_list)
    MCS_FS_Dicts[label]['Low']['Median Compressive Hist'], bin_edg,_ = scipy.stats.binned_statistic(all_fs_low['Frequency'], all_fs_low['Compressive Power'], statistic='median', bins=bin_list)
    
    MCS_FS_Dicts[label]['High']['Median Transverse Hist'], bin_edg,_ = scipy.stats.binned_statistic(all_fs_high['Frequency'], all_fs_high['Total Transverse Power'], statistic='median', bins=bin_list)
    MCS_FS_Dicts[label]['High']['Median Compressive Hist'], bin_edg,_ = scipy.stats.binned_statistic(all_fs_high['Frequency'], all_fs_high['Compressive Power'], statistic='median', bins=bin_list)

    #find the mean proton and alpha predicted frequencies
    MCS_FS_Dicts[label]['Low']['Mean Takahashi Frequency'] = np.mean(np.array(MCS_Dicts[label]['Low']['Takahashi Frequency']))
    MCS_FS_Dicts[label]['High']['Mean Takahashi Frequency'] = np.mean(np.array(MCS_Dicts[label]['High']['Takahashi Frequency']))
    MCS_FS_Dicts[label]['Low']['Median Takahashi Frequency'] = np.median(np.array(MCS_Dicts[label]['Low']['Takahashi Frequency']))
    MCS_FS_Dicts[label]['High']['Median Takahashi Frequency'] = np.median(np.array(MCS_Dicts[label]['High']['Takahashi Frequency']))

    MCS_FS_Dicts[label]['Low']['Mean Peak Transverse Frequency'] = midpoint_list[np.argmax(MCS_FS_Dicts[label]['Low']['Mean Transverse Hist'])]
    MCS_FS_Dicts[label]['High']['Mean Peak Transverse Frequency'] = midpoint_list[np.argmax(MCS_FS_Dicts[label]['High']['Mean Transverse Hist'])]
    MCS_FS_Dicts[label]['Low']['Mean Peak Compressive Frequency'] = midpoint_list[np.argmax(MCS_FS_Dicts[label]['Low']['Mean Compressive Hist'])]
    MCS_FS_Dicts[label]['High']['Mean Peak Compressive Frequency'] = midpoint_list[np.argmax(MCS_FS_Dicts[label]['High']['Mean Compressive Hist'])]
    MCS_FS_Dicts[label]['Low']['Median Peak Transverse Frequency'] = midpoint_list[np.argmax(MCS_FS_Dicts[label]['Low']['Median Transverse Hist'])]
    MCS_FS_Dicts[label]['High']['Median Peak Transverse Frequency'] = midpoint_list[np.argmax(MCS_FS_Dicts[label]['High']['Median Transverse Hist'])]
    MCS_FS_Dicts[label]['Low']['Median Peak Compressive Frequency'] = midpoint_list[np.argmax(MCS_FS_Dicts[label]['Low']['Median Compressive Hist'])]
    MCS_FS_Dicts[label]['High']['Median Peak Compressive Frequency'] = midpoint_list[np.argmax(MCS_FS_Dicts[label]['High']['Median Compressive Hist'])]
             
spectra_overlay_plot_comparison(MCS_Dicts['MCS 1']['Low']['Normed FS List'], MCS_Dicts['MCS 1']['High']['Normed FS List'], CS_name=r'Multi-case Study 1, $\theta_{Bx}$ 0-20°', trans_ave_1=MCS_FS_Dicts['MCS 1']['Low']['Mean Transverse Hist'], comp_ave_1=MCS_FS_Dicts['MCS 1']['Low']['Mean Compressive Hist'], trans_ave_2=MCS_FS_Dicts['MCS 1']['High']['Median Transverse Hist'], comp_ave_2=MCS_FS_Dicts['MCS 1']['High']['Median Compressive Hist'], midpoints=midpoint_list, tak_freq_1=MCS_FS_Dicts['MCS 1']['Low']['Median Takahashi Frequency'], tak_freq_2=MCS_FS_Dicts['MCS 1']['High']['Median Takahashi Frequency'])
#spectra_overlay_plot_comparison(low_0_20_fs_dfs_normed_A, high_0_20_fs_dfs_normed_A, CS_name='Case Study 1 Mean, 0-20°', trans_ave_1=mean_trans_hist_low_0_20_A, comp_ave_1=mean_comp_hist_low_0_20_A, trans_ave_2=mean_trans_hist_high_0_20_A, comp_ave_2=mean_comp_hist_high_0_20_A, midpoints=midpoint_list, tak_freq_1=mean_tak_0_20_A_low, tak_freq_2=mean_tak_0_20_A_high)

spectra_overlay_plot_comparison(MCS_Dicts['MCS 2A']['Low']['Normed FS List'], MCS_Dicts['MCS 2A']['High']['Normed FS List'], CS_name=r'Multi-case Study 2A, $\theta_{Bx}$ 20-40°', trans_ave_1=MCS_FS_Dicts['MCS 2A']['Low']['Median Transverse Hist'], comp_ave_1=MCS_FS_Dicts['MCS 2A']['Low']['Median Compressive Hist'], trans_ave_2=MCS_FS_Dicts['MCS 2A']['High']['Median Transverse Hist'], comp_ave_2=MCS_FS_Dicts['MCS 2A']['High']['Median Compressive Hist'], midpoints=midpoint_list, tak_freq_1=MCS_FS_Dicts['MCS 2A']['Low']['Median Takahashi Frequency'], tak_freq_2=MCS_FS_Dicts['MCS 2A']['High']['Median Takahashi Frequency'])
spectra_overlay_plot_comparison(MCS_Dicts['MCS 2B']['Low']['Normed FS List'], MCS_Dicts['MCS 2B']['High']['Normed FS List'], CS_name=r'Multi-case Study 2B, $\theta_{Bx}$ 20-40°', trans_ave_1=MCS_FS_Dicts['MCS 2B']['Low']['Median Transverse Hist'], comp_ave_1=MCS_FS_Dicts['MCS 2B']['Low']['Median Compressive Hist'], trans_ave_2=MCS_FS_Dicts['MCS 2B']['High']['Median Transverse Hist'], comp_ave_2=MCS_FS_Dicts['MCS 2B']['High']['Median Compressive Hist'], midpoints=midpoint_list, tak_freq_1=MCS_FS_Dicts['MCS 2B']['Low']['Median Takahashi Frequency'], tak_freq_2=MCS_FS_Dicts['MCS 2B']['High']['Median Takahashi Frequency'])
spectra_overlay_plot_comparison(MCS_Dicts['MCS 2C']['Low']['Normed FS List'], MCS_Dicts['MCS 2C']['High']['Normed FS List'], CS_name=r'Multi-case Study 2C, $\theta_{Bx}$ 20-40°', trans_ave_1=MCS_FS_Dicts['MCS 2C']['Low']['Median Transverse Hist'], comp_ave_1=MCS_FS_Dicts['MCS 2C']['Low']['Median Compressive Hist'], trans_ave_2=MCS_FS_Dicts['MCS 2C']['High']['Median Transverse Hist'], comp_ave_2=MCS_FS_Dicts['MCS 2C']['High']['Median Compressive Hist'], midpoints=midpoint_list, tak_freq_1=MCS_FS_Dicts['MCS 2C']['Low']['Median Takahashi Frequency'], tak_freq_2=MCS_FS_Dicts['MCS 2C']['High']['Median Takahashi Frequency'])

