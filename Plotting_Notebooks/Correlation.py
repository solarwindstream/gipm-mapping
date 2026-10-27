#Script to do linear correlation of two time series

from scipy import stats
import datetime as dt
import pandas as pd
import numpy as np

from gipm_plotting_modules import model_calcs, draw_background, set_limits
from CorrelationModules import minvar_analysis, minvar_rot, max_corr_slide_1_min_gse, max_corr_slide_1_min_MVA, ts_overlay

############################
############################

def slice_match_mva_4mins(int_start: pd.DateTime, inp_df1: pd.DataFrame, inp_df_2: pd.DataFrame, two_scs: list[str]):

    int_end = int_start + dt.timedelta(minutes=4)
    inp_df1_4min=inp_df1.loc[((inp_df1.index >= int_start) & (inp_df1.index < int_end))]
    #MVA transform for inp1
    v_1, w_1, inp_df1_4min = minvar_analysis(inp_df1_4min)
    print(v_1)
    #now apply same transform to inp2, using inp1 MVA frame

    temp_int_start = int_start - dt.timedelta(minutes=0.5)
    temp_int_end = int_start + dt.timedelta(minutes=4.5)

    temp_ts_short_2=inp_df_2.loc[((inp_df_2.index >= temp_int_start) & (inp_df_2.index < temp_int_end))]
    
    temp_ts_short_2=minvar_rot(temp_ts_short_2, v_1)
    #and then sliiiiiiiiiiiide
    mr_mag_t_2, mr_l_t_2, mr_m_t_2, mr_n_t_2 = max_corr_slide_1_min_MVA(temp_ts_short_2, int_start, inp_df1_4min)

    #match time for max corr of max transverse eigenvector (l)
    inp_df2_4min = temp_ts_short_2.loc[((temp_ts_short_2.index >= mr_l_t_2) & (temp_ts_short_2.index < (mr_l_t_2 + dt.timedelta(minutes=4))))]

    #find matching windows, shorten to exact same length and then plot
    diff = -1*abs(inp_df1_4min.shape[0] - inp_df2_4min.shape[0])
    
    if inp_df1_4min.shape[0] < inp_df2_4min.shape[0]:
        print("Df2 length: ", inp_df2_4min.shape[0], " longer than Df1 by", diff)
        inp_df2_4min = inp_df2_4min.iloc[:diff]
        print(inp_df1_4min.shape[0], inp_df2_4min.shape[0])
    
    if inp_df1_4min.shape[0] > inp_df2_4min.shape[0]:
        print("Df1 length: ",inp_df1_4min.shape[0], " longer than Df2 by", diff)
        inp_df1_4min = inp_df1_4min.iloc[:diff]
    print("New Df1 length: ",inp_df1_4min.shape[0], " Df2 length:", inp_df2_4min.shape[0])

    dt_strf = int_start.strftime("%y_%m_%d_%H%M")
    
    ts_overlay(inp_df1_4min, inp_df2_4min, f"MVA_{two_scs[0]}{two_scs[1]}_{dt_strf}_Overlay", two_scs, props=['B_mag', 'B_l', 'B_m', 'B_n'])

############################
############################

int_csv_files =['/Users/roseatkinson/Documents/Cluster_Int_CSVs/2017-03-07 21:18:00_C1.csv','/Users/roseatkinson/Documents/Cluster_Int_CSVs/2017-03-07 22:34:00_C3.csv', '/Users/roseatkinson/Documents/Cluster_Int_CSVs/2017-03-07 22:06:00_C4.csv']

ts_files = ['/Users/roseatkinson/Documents/Cluster_Time_Series/TS_2017-03-07 21:00:00_C1.csv','/Users/roseatkinson/Documents/Cluster_Time_Series/TS_2017-03-07 21:00:00_C3.csv','/Users/roseatkinson/Documents/Cluster_Time_Series/TS_2017-03-07 21:00:00_C4.csv']

df_dict = {}

scs=['C1','C3','C4']

for sc, int_file, ts_file in zip(scs, int_csv_files,ts_files):
    df_dict[sc]={}
    df = pd.read_csv(int_file,encoding='utf-8')
    df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
    df.set_index('datetime', inplace = True)
    df['sc']=sc
    df_dict[sc]['Interval Data']=df
    
    df = pd.read_csv(ts_file,encoding='utf-8')
    df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
    df.set_index('datetime', inplace = True)
    df_dict[sc]['Time Series']=df

start_int = pd.to_datetime("2017-03-08 00:50:00")

ts_df_C1=df_dict['C1']['Time Series']
ts_df_C3=df_dict['C3']['Time Series']

slice_match_mva_4mins(start_int, ts_df_C1, ts_df_C3, ['C1','C3'])

int_locs = {}

for sc in scs:
    int_locs[sc]={}
    int_locs[sc]['Data']=df_dict[sc]['Interval Data'].loc[[start_int]]
    print(sc)
    for axis in ['GSE X (mean)', 'GSE Y (mean)', 'GSE Z (mean)']:
        print(f"{axis}: {int_locs[sc]['Data'].at[start_int, axis]}")

#calculate distance between two spacecraft

r_13 = np.sqrt((int_locs['C1']['Data'].at[start_int, 'GSE X (mean)']-int_locs['C3']['Data'].at[start_int, 'GSE X (mean)'])**2+(int_locs['C1']['Data'].at[start_int, 'GSE Y (mean)']-int_locs['C3']['Data'].at[start_int, 'GSE Y (mean)'])**2+(int_locs['C1']['Data'].at[start_int, 'GSE Z (mean)']-int_locs['C3']['Data'].at[start_int, 'GSE Z (mean)'])**2)

print(f"Distance between spacecraft is {r_13} R_GSE")

#now repeat for 5Hz

five_hz_dict={}

for sc in scs:
    five_hz_dict[sc]={}
    path=f'/Users/roseatkinson/Documents/5Hz_Spectra/2017-03-08 00:00:00.100000_{sc}.csv'
    df = pd.read_csv(path,encoding='utf-8')
    df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
    df.set_index('datetime', inplace = True)
    five_hz_dict[sc]['Time Series']=df

slice_match_mva_4mins(start_int, five_hz_dict['C1']['Time Series'], five_hz_dict['C3']['Time Series'], ['C1','C3'])

#now repeat for Archer et al
five_hz_dict_a={}
scs=['C1', 'C2', 'C3','C4']
for sc in scs:
    five_hz_dict_a[sc]={}
    path=f'/Users/roseatkinson/Documents/5Hz_Spectra/2003-02-18 13:00:00.100000_{sc}.csv'
    df = pd.read_csv(path,encoding='utf-8')
    df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
    df.set_index('datetime', inplace = True)
    five_hz_dict_a[sc]['Time Series']=df

print('\033[1mArcher 2005 interval\033[0m')
start_int_a=pd.to_datetime("2003-02-18 13:45:00")
slice_match_mva_4mins(start_int_a, five_hz_dict_a['C1']['Time Series'], five_hz_dict_a['C2']['Time Series'], ['C1','C2'])
slice_match_mva_4mins(start_int_a, five_hz_dict_a['C1']['Time Series'], five_hz_dict_a['C3']['Time Series'], ['C1','C3'])
slice_match_mva_4mins(start_int_a, five_hz_dict_a['C1']['Time Series'], five_hz_dict_a['C4']['Time Series'], ['C1','C4'])

