#Script to do linear correlation of two time series

from scipy import stats
import datetime as dt
import pandas as pd
import numpy as np

from gipm_plotting_modules import model_calcs, draw_background, set_limits
from CorrelationModules import minvar_analysis, max_corr_slide_1_min_gse, max_corr_slide_1_min_MVA, ts_overlay

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

#int_start = datetime.strftime("%Y-%m-%d %H:%M:%S")

int_start = pd.to_datetime("2017-03-08 00:50:00")
int_end = int_start + dt.timedelta(minutes=4)

#find matching time series for C1
ts_df_C1=df_dict[sc]['Time Series']
ts_short_C1=ts_df_C1.loc[((ts_df_C1.index >= int_start) & (ts_df_C1.index < int_end))]

#MVA transform for C1
vrot_1, v_1, w_1 = minvar_analysis(t_s=ts_df_C1)
print(v_1, w_1)

#now apply same transform to C3

temp_int_start = int_start - dt.timedelta(minutes=0.5)
temp_int_end = int_start + dt.timedelta(minutes=4.5)

ts_df_C3=df_dict['C3']['Time Series']
temp_ts_short_C3=ts_df_C3.loc[((ts_df_C3.index >= temp_int_start) & (ts_df_C3.index < temp_int_end))]

# ts_short={}

# for sc in scs:
#     ts_df=df_dict[sc]['Time Series']
#     ts_short[sc] = ts_df.loc[((ts_df.index >= int_start) & (ts_df.index < int_end))]

#Bx_gse, By_gse, Bz_gse, B_mag
#for sc in scs: print (ts_short[sc].index[0])

#max_corr_slide_1_min_gse('C3', int_start, ts_short['C1'])
                
#'GSE X (mean)' 'SW Vx (mean)'

int_locs = {}

for sc in scs:
    int_locs[sc]={}
    int_locs[sc]['Data']=df_dict[sc]['Interval Data'].loc[[int_start]]
    print(sc)
    for axis in ['GSE X (mean)', 'GSE Y (mean)', 'GSE Z (mean)']:
        print(f"{axis}: {int_locs[sc]['Data'].at[int_start, axis]}")

#calculate distance between two spacecraft

r_13 = np.sqrt((int_locs['C1']['Data'].at[int_start, 'GSE X (mean)']-int_locs['C3']['Data'].at[int_start, 'GSE X (mean)'])**2+(int_locs['C1']['Data'].at[int_start, 'GSE Y (mean)']-int_locs['C3']['Data'].at[int_start, 'GSE Y (mean)'])**2+(int_locs['C1']['Data'].at[int_start, 'GSE Z (mean)']-int_locs['C3']['Data'].at[int_start, 'GSE Z (mean)'])**2)

print(f"Distance between spacecraft is {r_13} R_GSE")

#print(C1_0050_data, C3_0050_data, C4_0050_data)
# print(C1_0050_data.at[int_start, 'GSE X (mean)'], C3_0050_data.at[int_start, 'GSE X (mean)'], C4_0050_data.at[int_start, 'GSE X (mean)'])
#print(C3_0050_data.loc[['GSE X (mean)', 'GSE Y (mean)', 'GSE Z (mean)']])
#print(C4_0050_data.loc[['GSE X (mean)', 'GSE Y (mean)', 'GSE Z (mean)']])

#get GSE loc for 00:50 C3, C4 cases 
#determine if going backward in time or forward in time from C1 by looking (mostly) at GSE X values.
#calculate expected advection time and then shift in a +-30s window around it to find max correlation value


# loc_x = []
# loc_y = []
# loc_z = []

# for sc in scs:
#     loc_x.append(int_locs[sc]['Data'].at[int_start, 'GIPM X (OMNI mean)'])
#     loc_y.append(int_locs[sc]['Data'].at[int_start, 'GIPM Y (OMNI mean)'])
#     loc_z.append(int_locs[sc]['Data'].at[int_start, 'GIPM Z (OMNI mean)'])

# time_series_comp(df_dict, ['C1', 'C3'],fname='C1C3_2017_03_08_00_50_Time_Series')
# loc_plot(loc_x[0:2],loc_y[0:2],fname='C1C3_2017_03_08_00_50_Locations')   


