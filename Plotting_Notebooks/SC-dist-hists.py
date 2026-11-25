#Script to do linear correlation of two time series

from scipy import stats
import datetime as dt
import pandas as pd
import numpy as np
import glob

from gipm_plotting_modules import model_calcs, draw_background, set_limits
from CorrelationModules import minvar_analysis, minvar_rot, max_corr_slide_1_min_gse, max_corr_slide_1_min_MVA, ts_overlay

 ##load integrated Cluster CSVs

path = "/Users/roseatkinson/Documents/Cluster_Int_CSVs/**"

cl_df_list = {'C1': [], 'C2': [], 'C3': [], 'C4': []}

for path in glob.glob(path, recursive=True):
    if '.csv' in path:
        if 'C1' in path:
            df = pd.read_csv(path,encoding='utf-8')
            df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
            df.set_index('datetime', inplace = True)
            df['sc']=sc
            cl_df_list['C1'].append(df)
        if 'C2' in path:
            df = pd.read_csv(path,encoding='utf-8')
            df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
            df.set_index('datetime', inplace = True)
            df['sc']=sc
            cl_df_list['C2'].append(df)
        if 'C3' in path:
            df = pd.read_csv(path,encoding='utf-8')
            df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
            df.set_index('datetime', inplace = True)
            df['sc']=sc
            cl_df_list['C3'].append(df)
        if 'C4' in path:
            df = pd.read_csv(path,encoding='utf-8')
            df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
            df.set_index('datetime', inplace = True)
            df['sc']=sc
            cl_df_list['C4'].append(df)

cl_dfs = {}
scs = ['C1', 'C2', 'C3', 'C4']

for sc in scs:
    cl_dfs[scs]=pd.concat(cl_df_list[sc])

print("dfs loaded in")
print("C1 dfs len:", cl_dfs['C1'].shape)
print("C2 dfs len:", cl_dfs['C2'].shape)
print("C3 dfs len:", cl_dfs['C3'].shape)
print("C4 dfs len:", cl_dfs['C4'].shape)

#now for all intervals in C1, check if matching C2/C3/C4 interval and calculate GSE distance.

#C1-C2 
C1_int_start_list_C2m = []
C1_C2_matching_time = []
C1_C2_distance = []

for start_time in cl_dfs['C1'].index:
    if start_time in cl_dfs['C2'].index:
        c1_x = cl_dfs['C1'].at[start_time, 'GSE X (mean)']
        c1_y = cl_dfs['C1'].at[start_time, 'GSE Y (mean)']
        c1_z = cl_dfs['C1'].at[start_time, 'GSE Z (mean)']
        c2_x = cl_dfs['C2'].at[start_time, 'GSE X (mean)']
        c2_y = cl_dfs['C2'].at[start_time, 'GSE Y (mean)']
        c2_z = cl_dfs['C2'].at[start_time, 'GSE Z (mean)']    
        
        r_gse = np.sqrt((c1_x-c2_x)**2+(c1_y-c2_y)**2+(c1_z-c2_z)**2)
        
        C1_int_start_list_C2m.append(start_time)
        C1_C2_matching_time.append(start_time)
        C1_C2_distance.append(r_gse)
        
        
    

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
