import cdflib
import datetime as dt
import pandas as pd
import numpy as np
import math
import statistics
import os
from itertools import product

from CDF_ProcessingScripts.Cluster_CDF_conv import Cluster_cdf_conv

from CDF_ProcessingScripts.window_det_4min import window_det
from CDF_ProcessingScripts.omni_seg_centered import omni_seg, omni_seg_ap_ratio
from CDF_ProcessingScripts.gipm_transform_coeffs import gipm_transform_coeffs_mean, gipm_transform_coeffs_median
from CDF_ProcessingScripts.GIPM_loc_conv import gipm_loc_transform
from CDF_ProcessingScripts.new_xyz import new_xyz
from CDF_ProcessingScripts.FFT_Hann import FFT_Hann

#select year and spacecraft

slurm_id = os.environ("SLURM_ARRAY_TASK_ID")

years = range(2001, 2025)
spacecraft = range(1,5)
all_combos = list(product(years,spacecraft)) # list allows indexed access

cur_folder = f"/data/SPCS-HIETALA-Shocks/GIPM-MAPPING/Cluster_CDFs/CDFs_{all_combos[slurm_id][0]}/C{all_combos[slurm_id][1]}/Full_CDFs"
all_files = os.listdir(cur_folder)

#load OMNI NaNp CSV
omni_nanp_csv_path = '/data/home/apx059/OMNI_23yrs_NaNp.csv'
omni_nanp_df = pd.read_csv(omni_nanp_csv_path,encoding='utf-8')
omni_nanp_df['datetime'] = pd.to_datetime(omni_nanp_df['datetime'], format='mixed')
omni_nanp_df.set_index('datetime', inplace = True)

#load main yearly OMNI CSV
year_n = all_combos[slurm_id][0]
omni_csv_path = '/data/SPCS-HIETALA-Shocks/GIPM-MAPPING/New_OMNI_Raw_Files/CSVs/Raw_OMNI_' + year_n + '.csv'
om = pd.read_csv(omni_csv_path,encoding='utf-8')
om['datetime'] = pd.to_datetime(omni['datetime'], format='mixed')
om.set_index('datetime', inplace = True)