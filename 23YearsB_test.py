#23 year B values. not yet w/regard to OMNI. just appending B min/max/mean. quick test first.

import datetime as dt
import pandas as pd
import cdflib

from C1_Cluster_CDF_conv import C1_cdf_conv
from C2_Cluster_CDF_conv import C2_cdf_conv
from C3_Cluster_CDF_conv import C3_cdf_conv
from C4_Cluster_CDF_conv import C4_cdf_conv

year_n = '2002'
cdf_list_path = '/data/scratch/apx059/23_Years_Data/CDFs/CDFs_'
list_of_cdfs = cdf_list_path + year_n + '.csv'
#raw cdfs into dfs!
list_files = pd.read_csv(list_of_cdfs, header=None)
#list_files_T = list_files.T
#print('list_files_T:', list_files_T)
#list_files_mask = list_files.loc[((list_files.index >= batch_start) & (list_files.index < batch_end))]
list_files = list_files.rename(columns={0:'fnames'})
list_cdfs = list_files['fnames'].to_list()
year_path = '/data/scratch/apx059/23_Years_Data/CDFs/CDFs_' + year_n + '/'

df_list_c1 = []

for i in list_cdfs:
    if 'C1' in i and len(df_list_c1)<3:
        fpath = year_path + i
        print(i)
        df_c1 = C1_cdf_conv(fpath)
        a = df_c1.empty
        if not a:
            df_list_c1.append(df_c1)
    elif len(df_list_c1)>=3:
        break            
        
#now also read in list of C1 csvs
CSV_list_path = '/data/scratch/apx059/23_Years_Data/CSVs/C1/list_of_files.csv'
sc_path = '/data/scratch/apx059/23_Years_Data/CSVs/C1/'
list_CSV_files = pd.read_csv(CSV_list_path, header=None)
list_CSV_files = list_CSV_files.rename(columns={0:'fnames'})
list_csvs = list_CSV_files['fnames'].to_list()
print(list_csvs[1342:1345])

gipm_df_c1 = []

for i in list_csvs[1342:1345]:
    fpath = sc_path + i
    df_c1 = pd.read_csv(fpath)
    df_c1['datetime'] = pd.to_datetime(df_c1['datetime'])
    df_c1 = df_c1.set_index('datetime')
    gipm_df_c1.append(df_c1)
    
print(gipm_df_c1[0])

#this bit is for the new dfs 

list_expanded_dfs = []

time_window = dt.timedelta(seconds=120)

for i,j in zip(df_list_c1, gipm_df_c1):
    
    cl_min_list = []
    cl_max_list = []
    cl_mean_list = []
    times = []
    
    for m in j.index:
        start_time = m
        end_time = m + time_window
        mask = i.loc[(i.index >= start_time) & (i.index < end_time)]
        Cluster_list = mask['B_mag'].tolist()
        if Cluster_list:
            Cluster_min = min(Cluster_list)
            Cluster_max = max(Cluster_list)
            Cluster_mean = sum(Cluster_list)/len(Cluster_list)
            cl_min_list.append(Cluster_min)
            cl_mean_list.append(Cluster_mean)
            cl_max_list.append(Cluster_max)
            times.append(m)
    
    if cl_min_list:
        print('cl_min_list', cl_min_list[0:7])
        B_val_df = pd.DataFrame({'datetime': times,'B min': cl_min_list, 'B mean': cl_mean_list, 'B max': cl_max_list})
        B_val_df = B_val_df.set_index('datetime')
        new_cl_df = j.join([B_val_df])
        list_expanded_dfs.append(new_cl_df)

#save in testing folder:

testing_file_path = '/data/scratch/apx059/23_Years_Data/Testing/'

for i in list_expanded_dfs:
    firstwin = i.index[0]
    firstwin = str(firstwin)
    fpath = testing_file_path + firstwin + 'GIPMandB_C1.csv'
    i.to_csv(fpath)