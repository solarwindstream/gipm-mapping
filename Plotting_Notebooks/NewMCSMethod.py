import pandas as pd
import numpy as np
import datetime as dt
import glob

import matplotlib 
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
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

#iterate over x, y foreshock strips:
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
    df = cl_filtered_20_40.loc[(cl_filtered_20_40['GIPM X (OMNI mean)'] >= x) & (cl_filtered_20_40['GIPM X (OMNI mean)'] < x_end) & (cl_filtered_20_40['GIPM Y (OMNI mean)'] >= y_bounds[0]) & (cl_filtered_20_40['GIPM Y (OMNI mean)'] < y_bounds[1])]
    fs_df_list_20_40.append(df)

#Background start times and mean background conditions for all possible foreshock intervals
fs_dfs_0_20 = pd.concat(fs_df_list_0_20)
fs_dfs_20_40 = pd.concat(fs_df_list_20_40)

high_alpha=0.08
low_alpha=0.005

cl_intervals_foreshock = {'0_20':{},'20_40': {}}

cl_intervals_foreshock['0_20']['High'] = fs_dfs_0_20.loc[(fs_dfs_0_20['SW Na/Np (mean)'] > high_alpha)]
cl_intervals_foreshock['0_20']['Low'] = fs_dfs_0_20.loc[(fs_dfs_0_20['SW Na/Np (mean)'] < low_alpha)]
cl_intervals_foreshock['20_40']['High'] = fs_dfs_20_40.loc[(fs_dfs_20_40['SW Na/Np (mean)'] > high_alpha)]
cl_intervals_foreshock['20_40']['Low'] = fs_dfs_20_40.loc[(fs_dfs_20_40['SW Na/Np (mean)'] < low_alpha)]

#drop any intervals that do not have associated extreme F.S., and save those fourier spectra dfs
#in a dictionary with the interval start as the key.

errorfiles = []

path = '/Users/roseatkinson/Documents/Extreme_Spectra/'

FS_dict = {}
TS_dict = {}
actual_foreshock_intervals={}
for ca, alpha_dict in cl_intervals_foreshock.items():
    actual_foreshock_intervals[ca]={}
    FS_dict[ca]={}
    TS_dict[ca]={}
    for alpha, df in alpha_dict.items():
        actual_foreshock_intervals[ca][alpha]={}
        FS_dict[ca][alpha]={'C1':{}, 'C2': {}, 'C3':{}, 'C4': {}}
        TS_dict[ca][alpha]={'C1':{}, 'C2': {}, 'C3':{}, 'C4': {}}
        new_df_c1=pd.DataFrame(columns = df.columns, index=pd.to_datetime([]))
        new_df_c2=pd.DataFrame(columns = df.columns, index=pd.to_datetime([]))
        new_df_c3=pd.DataFrame(columns = df.columns, index=pd.to_datetime([]))
        new_df_c4=pd.DataFrame(columns = df.columns, index=pd.to_datetime([]))
        for datetime, row in df.iterrows():
            int_start = datetime.strftime("%Y-%m-%d %H:%M:%S")
            sc_no = row['sc']
            filename='FS_'+int_start+'_'+sc_no+'.csv'
            filepath = path+filename
            #transform units to mHz from Hz, & reverse diff from Ftak-Fpeak to Fpeak-Ftak:
            row['Takahashi Transverse Difference']= -1000*row['Takahashi Transverse Difference']
            row['Takahashi Compressive Difference']= -1000*row['Takahashi Compressive Difference'] 
            row['Peak Transverse Frequency'] = 1000*row['Peak Transverse Frequency']
            row['Peak Compressive Frequency'] = 1000*row['Peak Compressive Frequency']
            try:
                f_df = pd.read_csv(filepath,encoding='utf-8')
                #normalise power spectrum:
                omni_mean = row['IMF B (mean)']
                bmean_bomni = row['B mean/Bomni']
                cluster_b_mean = bmean_bomni*omni_mean
                f_df['Total Transverse Power'] = f_df['Total Transverse Power']/(cluster_b_mean**2)
                f_df['Compressive Power'] = f_df['Compressive Power']/(cluster_b_mean**2)
                win_end = datetime + dt.timedelta(minutes=4)
                ts_df = cl_ts[sc_no]
                ts_df_filt = ts_df.loc[((ts_df.index >= datetime) & (ts_df.index < win_end))]
                if sc_no=='C1':
                    new_df_c1.loc[datetime] = row
                    FS_dict[ca][alpha]['C1'][int_start]=f_df
                    TS_dict[ca][alpha]['C1'][int_start]=ts_df_filt
                if sc_no=='C2':
                    new_df_c2.loc[datetime] = row
                    FS_dict[ca][alpha]['C2'][int_start]=f_df
                    TS_dict[ca][alpha]['C2'][int_start]=ts_df_filt
                if sc_no=='C3':
                    new_df_c3.loc[datetime] = row
                    FS_dict[ca][alpha]['C3'][int_start]=f_df
                    TS_dict[ca][alpha]['C3'][int_start]=ts_df_filt
                if sc_no=='C4':
                    new_df_c4.loc[datetime] = row  
                    FS_dict[ca][alpha]['C4'][int_start]=f_df
                    TS_dict[ca][alpha]['C4'][int_start]=ts_df_filt
            except:
                errorfiles.append(filename)
                
        actual_foreshock_intervals[ca][alpha]['C1']=new_df_c1
        actual_foreshock_intervals[ca][alpha]['C2']=new_df_c2
        actual_foreshock_intervals[ca][alpha]['C3']=new_df_c3
        actual_foreshock_intervals[ca][alpha]['C4']=new_df_c4
        

#Pc4/3 ratio:
#Calculate from FS dfs, add results to integrated dfs
for ca, alpha_dict in actual_foreshock_intervals.items():
    for alpha, sc_dict in alpha_dict.items():
        #iterate over spacecraft too
        for sc, df in sc_dict.items():
            #establish lists that will be added as columns to new dataframe
            ratio_Pc4_Pc3_trans=[]
            ratio_Pc4_Pc3_comp=[]
            for datetime in df.index:
                int_start = datetime.strftime("%Y-%m-%d %H:%M:%S")
                FS=FS_dict[ca][alpha][sc][int_start]
                P_perp_low_f, P_perp_high_f, P_para_low_f, P_para_high_f = split_band_power_int(FS)
                ratio_trans = P_perp_low_f/P_perp_high_f
                ratio_comp = P_para_low_f/P_para_high_f
                ratio_Pc4_Pc3_trans.append(ratio_trans)
                ratio_Pc4_Pc3_comp.append(ratio_comp)
            df['Pc4/Pc3 Ratio, Transverse']=ratio_Pc4_Pc3_trans
            df['Pc4/Pc3 Ratio, Compressive']=ratio_Pc4_Pc3_comp

# print(actual_foreshock_intervals['0_20']['High']['C1'].index)
# for key, val in FS_dict['0_20']['High']['C1'].items():
#     print(key) 
#now select foreshock locations 1, 2A,B,C

#1

foreshock_intervals_list_1={}
FS_dict_1 = {}

for alpha, sc_dict in actual_foreshock_intervals['0_20'].items():
    foreshock_intervals_list_1[alpha]={}
    FS_dict_1[alpha]={}
    for sc, df in sc_dict.items():
        FS_dict_1[alpha][sc]={}
        x0,x1 = 8, 13
        y0, y1 = -13,-9
        #find all relevant intervals
        filtered_df=df.loc[(df['GIPM X (OMNI mean)'] >= x0) & (df['GIPM X (OMNI mean)'] < x1) & (df['GIPM Y (OMNI mean)'] >= y0) & (df['GIPM Y (OMNI mean)'] < y1)]
        if len(filtered_df) != 0:
            foreshock_intervals_list_1[alpha][sc]=filtered_df
            for datetime in filtered_df.index:
                dt_str=datetime.strftime("%Y-%m-%d %H:%M:%S")
                FS_dict_1[alpha][sc][dt_str]=FS_dict['0_20'][alpha][sc][dt_str]


#2A

foreshock_intervals_list_2A={}
FS_dict_2A = {}

for alpha, sc_dict in actual_foreshock_intervals['20_40'].items():
    foreshock_intervals_list_2A[alpha]={}
    FS_dict_2A[alpha]={}
    for sc, df in sc_dict.items():
        FS_dict_2A[alpha][sc]={}
        x0,x1 = 14, 17
        y0, y1 = -2,0
        #find all relevant intervals
        filtered_df=df.loc[(df['GIPM X (OMNI mean)'] >= x0) & (df['GIPM X (OMNI mean)'] < x1) & (df['GIPM Y (OMNI mean)'] >= y0) & (df['GIPM Y (OMNI mean)'] < y1)]
        if len(filtered_df) != 0:
            foreshock_intervals_list_2A[alpha][sc]=filtered_df
            for datetime in filtered_df.index:
                dt_str=datetime.strftime("%Y-%m-%d %H:%M:%S")
                FS_dict_2A[alpha][sc][dt_str]=FS_dict['20_40'][alpha][sc][dt_str]

#2B

foreshock_intervals_list_2B={}
FS_dict_2B = {}

for alpha, sc_dict in actual_foreshock_intervals['20_40'].items():
    foreshock_intervals_list_2B[alpha]={}
    FS_dict_2B[alpha]={}
    for sc, df in sc_dict.items():
        FS_dict_2B[alpha][sc]={}
        x0,x1 = 14, 17
        y0, y1 = -4,-2
        #find all relevant intervals
        filtered_df=df.loc[(df['GIPM X (OMNI mean)'] >= x0) & (df['GIPM X (OMNI mean)'] < x1) & (df['GIPM Y (OMNI mean)'] >= y0) & (df['GIPM Y (OMNI mean)'] < y1)]
        if len(filtered_df) != 0:
            foreshock_intervals_list_2B[alpha][sc]=filtered_df
            for datetime in filtered_df.index:
                dt_str=datetime.strftime("%Y-%m-%d %H:%M:%S")
                FS_dict_2B[alpha][sc][dt_str]=FS_dict['20_40'][alpha][sc][dt_str]

foreshock_intervals_list_2C={}
FS_dict_2C = {}

for alpha, sc_dict in actual_foreshock_intervals['20_40'].items():
    foreshock_intervals_list_2C[alpha]={}
    FS_dict_2C[alpha]={}
    for sc, df in sc_dict.items():
        FS_dict_2C[alpha][sc]={}
        x0,x1 = 8, 12
        y0, y1 = -12,-9
        #find all relevant intervals
        filtered_df=df.loc[(df['GIPM X (OMNI mean)'] >= x0) & (df['GIPM X (OMNI mean)'] < x1) & (df['GIPM Y (OMNI mean)'] >= y0) & (df['GIPM Y (OMNI mean)'] < y1)]
        if len(filtered_df) != 0:
            foreshock_intervals_list_2C[alpha][sc]=filtered_df
            for datetime in filtered_df.index:
                dt_str=datetime.strftime("%Y-%m-%d %H:%M:%S")
                FS_dict_2C[alpha][sc][dt_str]=FS_dict['20_40'][alpha][sc][dt_str]

print('MCS filtered')

#Now: produce mapping & spectral overlays
#concatenate all spectra in a given category to find mean & median
list_low_FS_MCS_1=[]

for sc, sc_dict in FS_dict_1['Low'].items():
    for dt, df in sc_dict.items():
        print('Sc: ', sc, 'Int: ', dt)
        list_low_FS_MCS_1.append(df)

print('Low MCS1 Intervals:', len(list_low_FS_MCS_1))

list_high_FS_MCS_1=[]

for sc, sc_dict in FS_dict_1['High'].items():
    for dt, df in sc_dict.items():
        print('Sc: ', sc, 'Int: ', dt)
        list_high_FS_MCS_1.append(df)
    
print('High MCS1 Intervals:', len(list_high_FS_MCS_1))

all_low_MCS_1= pd.concat(list_low_FS_MCS_1)
all_high_MCS_1 = pd.concat(list_high_FS_MCS_1)

bin_list = np.arange(2.083333335e-03, 11.2, 4.16666667e-03)
bin_list=[0]+bin_list
print(bin_list)
#low: mean and median
mean_trans_low_MCS_1, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_1['Frequency'], all_low_MCS_1['Total Transverse Power'], statistic='mean', bins=bin_list)
mean_comp_low_MCS_1, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_1['Frequency'], all_low_MCS_1['Compressive Power'], statistic='mean', bins=bin_list)
median_trans_low_MCS_1, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_1['Frequency'], all_low_MCS_1['Total Transverse Power'], statistic='median', bins=bin_list)
median_comp_low_MCS_1, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_1['Frequency'], all_low_MCS_1['Compressive Power'], statistic='median', bins=bin_list)

#high: mean and median
mean_trans_high_MCS_1, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_1['Frequency'], all_high_MCS_1['Total Transverse Power'], statistic='mean', bins=bin_list)
mean_comp_high_MCS_1, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_1['Frequency'], all_high_MCS_1['Compressive Power'], statistic='mean', bins=bin_list)

median_trans_high_MCS_1, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_1['Frequency'], all_high_MCS_1['Total Transverse Power'], statistic='median', bins=bin_list)
median_comp_high_MCS_1, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_1['Frequency'], all_high_MCS_1['Compressive Power'], statistic='median', bins=bin_list)

#find the mean proton and alpha predicted frequencies - MCS 1 low
tak_MCS_1_low = []
pc4_pc3_trans_MCS_1_low = []
pc4_pc3_comp_MCS_1_low = []
tak_trans_diff_MCS_1_low = []
tak_comp_diff_MCS_1_low = []
  
for sc, sc_dict in foreshock_intervals_list_1['Low'].items():
    tak_MCS_1_low.extend(sc_dict['Takahashi Frequency, Hz'].to_list())
    pc4_pc3_trans_MCS_1_low.extend(sc_dict['Pc4/Pc3 Ratio, Transverse'].to_list())
    pc4_pc3_comp_MCS_1_low.extend(sc_dict['Pc4/Pc3 Ratio, Compressive'].to_list())
    tak_trans_diff_MCS_1_low.extend(sc_dict['Takahashi Transverse Difference'].to_list())
    tak_comp_diff_MCS_1_low.extend(sc_dict['Takahashi Compressive Difference'].to_list())
     
mean_tak_MCS_1_low = np.mean(np.array(tak_MCS_1_low))
median_tak_MCS_1_low = np.median(np.array(tak_MCS_1_low))

#MCS 1 high
tak_MCS_1_high = []
pc4_pc3_trans_MCS_1_high = []
pc4_pc3_comp_MCS_1_high = []
tak_trans_diff_MCS_1_high = []
tak_comp_diff_MCS_1_high = []

for sc, sc_dict in foreshock_intervals_list_1['High'].items():
    tak_MCS_1_high.extend(sc_dict['Takahashi Frequency, Hz'].to_list())
    pc4_pc3_trans_MCS_1_high.extend(sc_dict['Pc4/Pc3 Ratio, Transverse'].to_list())
    pc4_pc3_comp_MCS_1_high.extend(sc_dict['Pc4/Pc3 Ratio, Compressive'].to_list())
    tak_trans_diff_MCS_1_high.extend(sc_dict['Takahashi Transverse Difference'].to_list())
    tak_comp_diff_MCS_1_high.extend(sc_dict['Takahashi Compressive Difference'].to_list())
  
mean_tak_MCS_1_high = np.mean(np.array(tak_MCS_1_high))
median_tak_MCS_1_high = np.median(np.array(tak_MCS_1_high))

midpoint_list = np.arange(4.16666667e-03, 11.2, 4.16666667e-03)
spectra_overlay_plot_comparison(list_low_FS_MCS_1, list_high_FS_MCS_1, CS_name='Case Study 1, 0-20°, Median', trans_ave_1=median_trans_low_MCS_1, comp_ave_1=median_comp_low_MCS_1, trans_ave_2=median_trans_high_MCS_1, comp_ave_2=median_comp_high_MCS_1, midpoints=midpoint_list, tak_freq_1=median_tak_MCS_1_low, tak_freq_2=median_tak_MCS_1_high)

spectra_overlay_plot_comparison(list_low_FS_MCS_1, list_high_FS_MCS_1, CS_name='Case Study 1, 0-20°, Mean', trans_ave_1=mean_trans_low_MCS_1, comp_ave_1=mean_comp_low_MCS_1, trans_ave_2=mean_trans_high_MCS_1, comp_ave_2=mean_comp_high_MCS_1, midpoints=midpoint_list, tak_freq_1=mean_tak_MCS_1_low, tak_freq_2=mean_tak_MCS_1_high)


#Produce mapping & spectral overlays for 2A
#concatenate all spectra in a given category to find mean & median
list_low_FS_MCS_2A=[]

for sc, sc_dict in FS_dict_2A['Low'].items():
    for dt, df in sc_dict.items():
        print('Sc: ', sc, 'Int: ', dt)
        list_low_FS_MCS_2A.append(df)

print('Low MCS 2A Intervals:', len(list_low_FS_MCS_2A))

list_high_FS_MCS_2A=[]

for sc, sc_dict in FS_dict_2A['High'].items():
    for dt, df in sc_dict.items():
        print('Sc: ', sc, 'Int: ', dt)
        list_high_FS_MCS_2A.append(df)
    
print('High MCS 2A Intervals:', len(list_high_FS_MCS_2A))

all_low_MCS_2A= pd.concat(list_low_FS_MCS_2A)
all_high_MCS_2A = pd.concat(list_high_FS_MCS_2A)

bin_list = np.arange(2.083333335e-03, 11.2, 4.16666667e-03)
bin_list=[0]+bin_list
print(bin_list)
#low: mean and median
mean_trans_low_MCS_2A, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_2A['Frequency'], all_low_MCS_2A['Total Transverse Power'], statistic='mean', bins=bin_list)
mean_comp_low_MCS_2A, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_2A['Frequency'], all_low_MCS_2A['Compressive Power'], statistic='mean', bins=bin_list)
median_trans_low_MCS_2A, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_2A['Frequency'], all_low_MCS_2A['Total Transverse Power'], statistic='median', bins=bin_list)
median_comp_low_MCS_2A, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_2A['Frequency'], all_low_MCS_2A['Compressive Power'], statistic='median', bins=bin_list)

#high: mean and median
mean_trans_high_MCS_2A, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_2A['Frequency'], all_high_MCS_2A['Total Transverse Power'], statistic='mean', bins=bin_list)
mean_comp_high_MCS_2A, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_2A['Frequency'], all_high_MCS_2A['Compressive Power'], statistic='mean', bins=bin_list)

median_trans_high_MCS_2A, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_2A['Frequency'], all_high_MCS_2A['Total Transverse Power'], statistic='median', bins=bin_list)
median_comp_high_MCS_2A, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_2A['Frequency'], all_high_MCS_2A['Compressive Power'], statistic='median', bins=bin_list)

#find the mean proton and alpha predicted frequencies - MCS 1 low
tak_MCS_2A_low = []
pc4_pc3_trans_MCS_2A_low = []
pc4_pc3_comp_MCS_2A_low = []
tak_trans_diff_MCS_2A_low = []
tak_comp_diff_MCS_2A_low = []
peak_f_t_MCS_2A_low = []
peak_f_c_MCS_2A_low = []

for sc, sc_dict in foreshock_intervals_list_2A['Low'].items():
    tak_MCS_2A_low.extend(sc_dict['Takahashi Frequency, Hz'].to_list())
    pc4_pc3_trans_MCS_2A_low.extend(sc_dict['Pc4/Pc3 Ratio, Transverse'].to_list())
    pc4_pc3_comp_MCS_2A_low.extend(sc_dict['Pc4/Pc3 Ratio, Compressive'].to_list())
    tak_trans_diff_MCS_2A_low.extend(sc_dict['Takahashi Transverse Difference'].to_list())
    tak_comp_diff_MCS_2A_low.extend(sc_dict['Takahashi Compressive Difference'].to_list())
    peak_f_t_MCS_2A_low.extend(sc_dict['Peak Transverse Frequency'].to_list())
    peak_f_c_MCS_2A_low.extend(sc_dict['Peak Compressive Frequency'].to_list())
    
mean_tak_MCS_2A_low = np.mean(np.array(tak_MCS_2A_low))
median_tak_MCS_2A_low = np.median(np.array(tak_MCS_2A_low))

#MCS 2A high
tak_MCS_2A_high = []
pc4_pc3_trans_MCS_2A_high = []
pc4_pc3_comp_MCS_2A_high = []
tak_trans_diff_MCS_2A_high = []
tak_comp_diff_MCS_2A_high = []
peak_f_t_MCS_2A_high = []
peak_f_c_MCS_2A_high = []

for sc, sc_dict in foreshock_intervals_list_2A['High'].items():
    tak_MCS_2A_high.extend(sc_dict['Takahashi Frequency, Hz'].to_list())
    pc4_pc3_trans_MCS_2A_high.extend(sc_dict['Pc4/Pc3 Ratio, Transverse'].to_list())
    pc4_pc3_comp_MCS_2A_high.extend(sc_dict['Pc4/Pc3 Ratio, Compressive'].to_list())
    tak_trans_diff_MCS_2A_high.extend(sc_dict['Takahashi Transverse Difference'].to_list())
    tak_comp_diff_MCS_2A_high.extend(sc_dict['Takahashi Compressive Difference'].to_list())
    peak_f_t_MCS_2A_high.extend(sc_dict['Peak Transverse Frequency'].to_list())
    peak_f_c_MCS_2A_high.extend(sc_dict['Peak Compressive Frequency'].to_list())
        
mean_tak_MCS_2A_high = np.mean(np.array(tak_MCS_2A_high))
median_tak_MCS_2A_high = np.median(np.array(tak_MCS_2A_high))

midpoint_list = np.arange(4.16666667e-03, 11.2, 4.16666667e-03)
spectra_overlay_plot_comparison(list_low_FS_MCS_2A, list_high_FS_MCS_2A, CS_name='Case Study 2A, 20-40°, Median', trans_ave_1=median_trans_low_MCS_2A, comp_ave_1=median_comp_low_MCS_2A, trans_ave_2=median_trans_high_MCS_2A, comp_ave_2=median_comp_high_MCS_2A, midpoints=midpoint_list, tak_freq_1=median_tak_MCS_2A_low, tak_freq_2=median_tak_MCS_2A_high)

spectra_overlay_plot_comparison(list_low_FS_MCS_2A, list_high_FS_MCS_2A, CS_name='Case Study 2A, 20-40°, Mean', trans_ave_1=mean_trans_low_MCS_2A, comp_ave_1=mean_comp_low_MCS_2A, trans_ave_2=mean_trans_high_MCS_2A, comp_ave_2=mean_comp_high_MCS_2A, midpoints=midpoint_list, tak_freq_1=mean_tak_MCS_2A_low, tak_freq_2=mean_tak_MCS_2A_high)

#Now: produce mapping & spectral overlays for 2B
#concatenate all spectra in a given category to find mean & mean
list_low_FS_MCS_2B=[]

for sc, sc_dict in FS_dict_2B['Low'].items():
    for dt, df in sc_dict.items():
        print('Sc: ', sc, 'Int: ', dt)
        list_low_FS_MCS_2B.append(df)

print('Low MCS1 Intervals:', len(list_low_FS_MCS_2B))

list_high_FS_MCS_2B=[]

for sc, sc_dict in FS_dict_2B['High'].items():
    for dt, df in sc_dict.items():
        print('Sc: ', sc, 'Int: ', dt)
        list_high_FS_MCS_2B.append(df)
    
print('High MCS1 Intervals:', len(list_high_FS_MCS_2B))

all_low_MCS_2B= pd.concat(list_low_FS_MCS_2B)
all_high_MCS_2B = pd.concat(list_high_FS_MCS_2B)

bin_list = np.arange(2.083333335e-03, 11.2, 4.16666667e-03)
bin_list=[0]+bin_list
print(bin_list)
#low: mean and median
mean_trans_low_MCS_2B, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_2B['Frequency'], all_low_MCS_2B['Total Transverse Power'], statistic='mean', bins=bin_list)
mean_comp_low_MCS_2B, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_2B['Frequency'], all_low_MCS_2B['Compressive Power'], statistic='mean', bins=bin_list)
median_trans_low_MCS_2B, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_2B['Frequency'], all_low_MCS_2B['Total Transverse Power'], statistic='median', bins=bin_list)
median_comp_low_MCS_2B, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_2B['Frequency'], all_low_MCS_2B['Compressive Power'], statistic='median', bins=bin_list)

#high: mean and median
mean_trans_high_MCS_2B, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_2B['Frequency'], all_high_MCS_2B['Total Transverse Power'], statistic='mean', bins=bin_list)
mean_comp_high_MCS_2B, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_2B['Frequency'], all_high_MCS_2B['Compressive Power'], statistic='mean', bins=bin_list)

median_trans_high_MCS_2B, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_2B['Frequency'], all_high_MCS_2B['Total Transverse Power'], statistic='median', bins=bin_list)
median_comp_high_MCS_2B, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_2B['Frequency'], all_high_MCS_2B['Compressive Power'], statistic='median', bins=bin_list)

#find the mean proton and alpha predicted frequencies - MCS 2B low
tak_MCS_2B_low = []
pc4_pc3_trans_MCS_2B_low = []
pc4_pc3_comp_MCS_2B_low = []
tak_trans_diff_MCS_2B_low = []
tak_comp_diff_MCS_2B_low = []
peak_f_t_MCS_2B_low = []
peak_f_c_MCS_2B_low = []
    
for sc, sc_dict in foreshock_intervals_list_2B['Low'].items():
    tak_MCS_2B_low.extend(sc_dict['Takahashi Frequency, Hz'].to_list())
    pc4_pc3_trans_MCS_2B_low.extend(sc_dict['Pc4/Pc3 Ratio, Transverse'].to_list())
    pc4_pc3_comp_MCS_2B_low.extend(sc_dict['Pc4/Pc3 Ratio, Compressive'].to_list())
    tak_trans_diff_MCS_2B_low.extend(sc_dict['Takahashi Transverse Difference'].to_list())
    tak_comp_diff_MCS_2B_low.extend(sc_dict['Takahashi Compressive Difference'].to_list())
    peak_f_t_MCS_2B_low.extend(sc_dict['Peak Transverse Frequency'].to_list())
    peak_f_c_MCS_2B_low.extend(sc_dict['Peak Compressive Frequency'].to_list())
    
mean_tak_MCS_2B_low = np.mean(np.array(tak_MCS_2B_low))
median_tak_MCS_2B_low = np.median(np.array(tak_MCS_2B_low))

#MCS 2B high
tak_MCS_2B_high = []
pc4_pc3_trans_MCS_2B_high = []
pc4_pc3_comp_MCS_2B_high = []
tak_trans_diff_MCS_2B_high = []
tak_comp_diff_MCS_2B_high = []
peak_f_t_MCS_2B_high = []
peak_f_c_MCS_2B_high = []

for sc, sc_dict in foreshock_intervals_list_2B['High'].items():
    tak_MCS_2B_high.extend(sc_dict['Takahashi Frequency, Hz'].to_list())
    pc4_pc3_trans_MCS_2B_high.extend(sc_dict['Pc4/Pc3 Ratio, Transverse'].to_list())
    pc4_pc3_comp_MCS_2B_high.extend(sc_dict['Pc4/Pc3 Ratio, Compressive'].to_list())
    tak_trans_diff_MCS_2B_high.extend(sc_dict['Takahashi Transverse Difference'].to_list())
    tak_comp_diff_MCS_2B_high.extend(sc_dict['Takahashi Compressive Difference'].to_list())
    peak_f_t_MCS_2B_high.extend(sc_dict['Peak Transverse Frequency'].to_list())
    peak_f_c_MCS_2B_high.extend(sc_dict['Peak Compressive Frequency'].to_list())
    
mean_tak_MCS_2B_high = np.mean(np.array(tak_MCS_2B_high))
median_tak_MCS_2B_high = np.median(np.array(tak_MCS_2B_high))

midpoint_list = np.arange(4.16666667e-03, 11.2, 4.16666667e-03)
spectra_overlay_plot_comparison(list_low_FS_MCS_2B, list_high_FS_MCS_2B, CS_name='Case Study 2B, 20-40°, Median', trans_ave_1=median_trans_low_MCS_2B, comp_ave_1=median_comp_low_MCS_2B, trans_ave_2=median_trans_high_MCS_2B, comp_ave_2=median_comp_high_MCS_2B, midpoints=midpoint_list, tak_freq_1=median_tak_MCS_2B_low, tak_freq_2=median_tak_MCS_2B_high)

spectra_overlay_plot_comparison(list_low_FS_MCS_2B, list_high_FS_MCS_2B, CS_name='Case Study 2B, 20-40°, Mean', trans_ave_1=mean_trans_low_MCS_2B, comp_ave_1=mean_comp_low_MCS_2B, trans_ave_2=mean_trans_high_MCS_2B, comp_ave_2=mean_comp_high_MCS_2B, midpoints=midpoint_list, tak_freq_1=mean_tak_MCS_2B_low, tak_freq_2=mean_tak_MCS_2B_high)

#Now: produce mapping & spectral overlays
#concatenate all spectra in a given category to find mean & median
list_low_FS_MCS_2C=[]

for sc, sc_dict in FS_dict_2C['Low'].items():
    for dt, df in sc_dict.items():
        print('Sc: ', sc, 'Int: ', dt)
        list_low_FS_MCS_2C.append(df)

print('Low MCS1 Intervals:', len(list_low_FS_MCS_2C))

list_high_FS_MCS_2C=[]

for sc, sc_dict in FS_dict_2C['High'].items():
    for dt, df in sc_dict.items():
        print('Sc: ', sc, 'Int: ', dt)
        list_high_FS_MCS_2C.append(df)
    
print('High MCS1 Intervals:', len(list_high_FS_MCS_2C))

all_low_MCS_2C= pd.concat(list_low_FS_MCS_2C)
all_high_MCS_2C = pd.concat(list_high_FS_MCS_2C)

bin_list = np.arange(2.083333335e-03, 11.2, 4.16666667e-03)
bin_list=[0]+bin_list
print(bin_list)
#low: mean and median
mean_trans_low_MCS_2C, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_2C['Frequency'], all_low_MCS_2C['Total Transverse Power'], statistic='mean', bins=bin_list)
mean_comp_low_MCS_2C, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_2C['Frequency'], all_low_MCS_2C['Compressive Power'], statistic='mean', bins=bin_list)
median_trans_low_MCS_2C, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_2C['Frequency'], all_low_MCS_2C['Total Transverse Power'], statistic='median', bins=bin_list)
median_comp_low_MCS_2C, bin_edg,_ = scipy.stats.binned_statistic(all_low_MCS_2C['Frequency'], all_low_MCS_2C['Compressive Power'], statistic='median', bins=bin_list)

#high: mean and median
mean_trans_high_MCS_2C, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_2C['Frequency'], all_high_MCS_2C['Total Transverse Power'], statistic='mean', bins=bin_list)
mean_comp_high_MCS_2C, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_2C['Frequency'], all_high_MCS_2C['Compressive Power'], statistic='mean', bins=bin_list)

median_trans_high_MCS_2C, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_2C['Frequency'], all_high_MCS_2C['Total Transverse Power'], statistic='median', bins=bin_list)
median_comp_high_MCS_2C, bin_edg,_ = scipy.stats.binned_statistic(all_high_MCS_2C['Frequency'], all_high_MCS_2C['Compressive Power'], statistic='median', bins=bin_list)

#find the mean proton and alpha predicted frequencies - MCS 1 low
tak_MCS_2C_low = []
pc4_pc3_trans_MCS_2C_low = []
pc4_pc3_comp_MCS_2C_low = []
tak_trans_diff_MCS_2C_low = []
tak_comp_diff_MCS_2C_low=[]
peak_f_t_MCS_2C_low = []
peak_f_c_MCS_2C_low = []

for sc, sc_dict in foreshock_intervals_list_2C['Low'].items():
    tak_MCS_2C_low.extend(sc_dict['Takahashi Frequency, Hz'].to_list())
    pc4_pc3_trans_MCS_2C_low.extend(sc_dict['Pc4/Pc3 Ratio, Transverse'].to_list())
    pc4_pc3_comp_MCS_2C_low.extend(sc_dict['Pc4/Pc3 Ratio, Compressive'].to_list())
    tak_trans_diff_MCS_2C_low.extend(sc_dict['Takahashi Transverse Difference'].to_list())
    tak_comp_diff_MCS_2C_low.extend(sc_dict['Takahashi Compressive Difference'].to_list())
    peak_f_t_MCS_2C_low.extend(sc_dict['Peak Transverse Frequency'].to_list())
    peak_f_c_MCS_2C_low.extend(sc_dict['Peak Compressive Frequency'].to_list())
    
mean_tak_MCS_2C_low = np.mean(np.array(tak_MCS_2C_low))
median_tak_MCS_2C_low = np.median(np.array(tak_MCS_2C_low))

#MCS 2C high
tak_MCS_2C_high = []
pc4_pc3_trans_MCS_2C_high = []
pc4_pc3_comp_MCS_2C_high = []
tak_trans_diff_MCS_2C_high = []
tak_comp_diff_MCS_2C_high = []
peak_f_t_MCS_2C_high = []
peak_f_c_MCS_2C_high = []
    
for sc, sc_dict in foreshock_intervals_list_2C['High'].items():
    tak_MCS_2C_high.extend(sc_dict['Takahashi Frequency, Hz'].to_list())
    pc4_pc3_trans_MCS_2C_high.extend(sc_dict['Pc4/Pc3 Ratio, Transverse'].to_list())
    pc4_pc3_comp_MCS_2C_high.extend(sc_dict['Pc4/Pc3 Ratio, Compressive'].to_list())
    tak_trans_diff_MCS_2C_high.extend(sc_dict['Takahashi Transverse Difference'].to_list())
    tak_comp_diff_MCS_2C_high.extend(sc_dict['Takahashi Compressive Difference'].to_list())
    peak_f_t_MCS_2C_high.extend(sc_dict['Peak Transverse Frequency'].to_list())
    peak_f_c_MCS_2C_high.extend(sc_dict['Peak Compressive Frequency'].to_list())
        
mean_tak_MCS_2C_high = np.mean(np.array(tak_MCS_2C_high))
median_tak_MCS_2C_high = np.median(np.array(tak_MCS_2C_high))

midpoint_list = np.arange(4.16666667e-03, 11.2, 4.16666667e-03)
spectra_overlay_plot_comparison(list_low_FS_MCS_2C, list_high_FS_MCS_2C, CS_name='Case Study 2C, 20-40°, Median', trans_ave_1=median_trans_low_MCS_2C, comp_ave_1=median_comp_low_MCS_2C, trans_ave_2=median_trans_high_MCS_2C, comp_ave_2=median_comp_high_MCS_2C, midpoints=midpoint_list, tak_freq_1=median_tak_MCS_2C_low, tak_freq_2=median_tak_MCS_2C_high)

spectra_overlay_plot_comparison(list_low_FS_MCS_2C, list_high_FS_MCS_2C, CS_name='Case Study 2C, 20-40°, Mean', trans_ave_1=mean_trans_low_MCS_2C, comp_ave_1=mean_comp_low_MCS_2C, trans_ave_2=mean_trans_high_MCS_2C, comp_ave_2=mean_comp_high_MCS_2C, midpoints=midpoint_list, tak_freq_1=mean_tak_MCS_2C_low, tak_freq_2=mean_tak_MCS_2C_high)

################################
################################

#produce plot showing histograms for fpeak-ftak

#first produce for joint low + high 
# tak_trans_diff_MCS_1_high
# tak_comp_diff_MCS_1_high
low_tak_diff_trans = tak_trans_diff_MCS_1_low+tak_trans_diff_MCS_2A_low+tak_trans_diff_MCS_2B_low+tak_trans_diff_MCS_2C_low
low_tak_diff_comp = tak_comp_diff_MCS_1_low+tak_comp_diff_MCS_2A_low+tak_comp_diff_MCS_2B_low+tak_comp_diff_MCS_2C_low
high_tak_diff_trans = tak_trans_diff_MCS_1_high+tak_trans_diff_MCS_2A_high+tak_trans_diff_MCS_2B_high+tak_trans_diff_MCS_2C_high
high_tak_diff_comp = tak_comp_diff_MCS_1_high+tak_comp_diff_MCS_2A_high+tak_comp_diff_MCS_2B_high+tak_comp_diff_MCS_2C_high

low_mean_t_diff = np.mean(np.array(low_tak_diff_trans))
high_mean_t_diff = np.mean(np.array(high_tak_diff_trans))
# print(f"{high_tak_diff_trans=}")
# print(f"{high_mean_t_diff=}")
low_mean_c_diff= np.mean(np.array(low_tak_diff_comp))
high_mean_c_diff = np.mean(np.array(high_tak_diff_comp))
# print(f"{high_tak_diff_comp=}")
# print(f"{high_mean_c_diff=}")

bins=np.array([-35.4166667 , -31.25000003, -27.08333336, -22.91666669,
       -18.75000002, -14.58333335, -10.41666668,  -6.25000001,
        -2.08333334, 2.08333334,  6.25000001, 10.41666668, 14.58333335,
       18.75000002, 22.91666669, 27.08333336, 31.25000003,
       35.4166667])
# bins = 1000*np.array([-0.01458333, -0.01041667, -0.00625, -0.00208333, 0.00208333, 0.00625,0.01041667, 0.01458333])

fig = plt.figure(figsize=(9, 5))
gs = GridSpec(1, 2, figure=fig)

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1], sharey=ax1)
plt.setp(ax2.get_yticklabels(), visible=False)

ax1.hist(low_tak_diff_trans, bins=bins, histtype='barstacked', color='blue', ec='darkblue', linewidth=1, alpha=0.5, label='Low Alpha')
ax1.hist(high_tak_diff_trans, bins=bins, histtype='barstacked', color='orange', ec='darkorange', linewidth=1, alpha=0.5, label='High Alpha')
ax1.vlines(low_mean_t_diff, 0, 10, color='k', linestyles='-', label='Mean (Low Alpha)')
ax1.vlines(high_mean_t_diff, 0, 10, color='k', linestyles='dotted', label='Mean (High Alpha)')
#ax1.legend()
ax1.set_xlabel('Fpeak - Ftak, mHz')
ax1.set_ylabel('No. Intervals')
ax1.set_title('Transverse', fontsize=16)
ax1.set_xlim(-30,30)
ax1.set_ylim(0,10)

ax2.hist(low_tak_diff_comp, bins=bins, histtype='barstacked', color='blue', ec='darkblue', linewidth=1,  alpha=0.5, label='Low Alpha')
ax2.hist(high_tak_diff_comp, bins=bins, histtype='barstacked', color='orange', ec='darkorange', linewidth=1, alpha=0.5, label='High Alpha')
ax2.set_xlabel('Fpeak - Ftak, mHz')
ax2.set_title('Compressive', fontsize=16)
ax2.vlines(low_mean_c_diff, 0, 10, color='k', linestyles='-', label='Mean (Low Alpha)')
ax2.vlines(high_mean_c_diff, 0, 10, color='k', linestyles='dotted', label='Mean (High Alpha)')
ax2.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
ax2.set_xlim(-30,30)
ax2.set_ylim(0,10)

path = "/Users/roseatkinson/Documents/New_Figs/Fpeak_Ftak_New.png"
# if 'filename' in kwargs:
#     path = f"/Users/roseatkinson/Documents/New_Figs/{kwargs['filename']}.png"
# else:
#     path = "/Users/roseatkinson/Documents/New_Figs/AlphaCaseStudyLocations.png"
plt.savefig(path, bbox_inches='tight')

############
#Scatter Plot

list_low_tak = []

for l in [tak_MCS_2A_low,tak_MCS_2B_low, tak_MCS_2C_low]:
    ar = np.array(l)
    ar=ar*1000
    list_low_tak.append(ar)

list_high_tak = []

for l in [tak_MCS_2A_high,tak_MCS_2B_high, tak_MCS_2C_high]:
    ar = np.array(l)
    ar=ar*1000
    list_high_tak.append(ar) 
 
list_low_t = [peak_f_t_MCS_2A_low, peak_f_t_MCS_2B_low, peak_f_t_MCS_2C_low]

list_high_t = [peak_f_t_MCS_2A_high, peak_f_t_MCS_2B_high, peak_f_t_MCS_2C_high]
       
list_low_c = [peak_f_c_MCS_2A_low, peak_f_c_MCS_2B_low, peak_f_c_MCS_2C_low]

list_high_c = [peak_f_c_MCS_2A_high, peak_f_c_MCS_2B_high, peak_f_c_MCS_2C_high]
    
fig = plt.figure(figsize=(9, 5))
gs = GridSpec(1, 2, figure=fig)

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1], sharey=ax1)
plt.setp(ax2.get_yticklabels(), visible=False)

line=range(15,37,1)

#pick colors from matplotlib maps

# col_high = plt.cm.hot([0.8,0.6,0.4,0.2])
# col_low = plt.cm.winter_r([0,0.4,0.7,1])  

col_high = plt.cm.RdYlBu([0.4,0.2,0])
col_low = plt.cm.RdYlBu([0.6,0.8,1]) 

#ax1.legend()
ax1.plot(line, line, color='black')

for x,y,color in zip(list_low_tak, list_low_t, col_low):
    ax1.scatter(x, y,color=color, s=60)

for x,y,color in zip(list_high_tak, list_high_t, col_high):
    ax1.scatter(x, y,color=color, s=60)
    
ax1.set_xlabel('Ftak, mHz')
ax1.set_ylabel('Fpeak, mHz')
ax1.set_title('Transverse', fontsize=16)

ax2.plot(line, line, color='black')

labs_low = ['Low (2A)','Low (2B)', 'Low (2C)']
labs_high = ['High (2A)','High (2B)', 'High (2C)']

for x,y,color,lab in zip(list_low_tak, list_low_c, col_low, labs_low):
    ax2.scatter(x, y,color=color, s=60, label=lab)

for x,y,color,lab in zip(list_high_tak, list_high_c, col_high, labs_high):
    ax2.scatter(x, y,color=color, s=60, label=lab)
    
ax2.set_xlabel('Ftak, mHz')
ax2.set_title('Compressive', fontsize=16)

ax2.legend(bbox_to_anchor=(1.04, 1), loc="upper left")

path = "/Users/roseatkinson/Documents/New_Figs/Fpeak_Ftak_Scatter_New.png"
# if 'filename' in kwargs:
#     path = f"/Users/roseatkinson/Documents/New_Figs/{kwargs['filename']}.png"
# else:
#     path = "/Users/roseatkinson/Documents/New_Figs/AlphaCaseStudyLocations.png"
plt.savefig(path, bbox_inches='tight')