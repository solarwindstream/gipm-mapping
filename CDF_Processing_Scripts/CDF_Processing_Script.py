import cdflib
import datetime as dt
import pandas as pd
import numpy as np
import math
import os
from itertools import product

from cdfconv_gipm import cdfconv_gipm
from Cluster_CDF_conv import Cluster_cdf_conv

print('loaded modules', flush=True)
#select year and spacecraft

slurm_id = int(os.environ["SLURM_ARRAY_TASK_ID"])

years = range(2001, 2025)
spacecraft = range(1,5)
all_combos = list(product(years,spacecraft)) # list allows indexed access

cur_folder = f"/data/SPCS-HIETALA-Shocks/GIPM-MAPPING/Cluster_CDFs/CDFs_{all_combos[slurm_id][0]}/Full_CDFs/C{all_combos[slurm_id][1]}"
all_files = os.listdir(cur_folder)

print(f'Analysing folder {cur_folder}', flush=True)

#load OMNI NaNp CSV
omni_nanp_csv_path = '/data/home/apx059/OMNI_23yrs_NaNp.csv'
omni_nanp_df = pd.read_csv(omni_nanp_csv_path,encoding='utf-8')
omni_nanp_df['datetime'] = pd.to_datetime(omni_nanp_df['datetime'], format='mixed')
omni_nanp_df.set_index('datetime', inplace = True)
print('loaded OMNI NaNp', flush=True)

#load main yearly OMNI CSV, taking year from job ID.
year_n = str(all_combos[slurm_id][0])
omni_csv_path = '/data/SPCS-HIETALA-Shocks/GIPM-MAPPING/New_OMNI_Raw_Files/CSVs/Raw_OMNI_' + year_n + '.csv'
om = pd.read_csv(omni_csv_path,encoding='utf-8')
om['datetime'] = pd.to_datetime(om['datetime'], format='mixed')
om.set_index('datetime', inplace = True)

print('Loaded in OMNI raw', flush=True)

#load in cdfs & create list of dataframes of raw Cluster data from CDF file list!
sc_ref = 'C' + str(all_combos[slurm_id][1])
path = cur_folder + '/'
    
df_list = []
err_list: list[str] = []

for count, file in enumerate(all_files):
    cdf_path = path + file
    df = Cluster_cdf_conv(cdf_path, sc_ref)
    if not df.empty:
        df_list.append(df)
    else:
        err_list.append(file)
    print(count+1, 'out of', len(all_files), 'CDFs loaded in', flush=True)

print(f'Loaded in CDFs ({len(df_list)}/{len(all_files)} successful)', flush=True)
print('Error files:', flush=True)
for ef in err_list:
    print(ef, flush=True)
#input raw Cluster dfs into GIPM conversion module!

cdfconv_gipm(year_n, df_list, sc_ref, om, omni_nanp_df)


print('love you sasha <3', flush=True)
