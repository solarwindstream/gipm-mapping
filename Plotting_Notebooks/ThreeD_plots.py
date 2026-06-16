#Split location-only results by MA and cone angle, then go into those bins and pull out some spectra

import pandas as pd
import numpy as np
import datetime as dt
import glob
from heatmap_modules import compute_hists2d_low_data

from heatmap_plotting_modules import CA_MA_Binned_Plot
from Three_D_Modules import  ThreeD_Binned_Plot

##load Cluster CSVs

cl_file_list = []

path = "/Users/roseatkinson/Documents/Cluster_Int_CSVs/**"

for path in glob.glob(path, recursive=True):
    if '.csv' in path:
        cl_file_list.append(path)
    
cl_dfs = []

for file in cl_file_list:
    df = pd.read_csv(file,encoding='utf-8')
    df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
    df.set_index('datetime', inplace = True)
    cl_dfs.append(df)

cl_power_all = pd.concat(cl_dfs)

print("dfs loaded in")

cl_filtered = cl_power_all.loc[(cl_power_all['OMNI Dist from X line (mean)'] < 70) & (cl_power_all['Max IMF Deviation'] < 60)]
cl_filtered_0 = cl_filtered.loc[(cl_filtered['GIPM Z (OMNI mean)'] < 1) & (cl_filtered['GIPM Z (OMNI mean)'] > -1)]
cl_filtered_1 = cl_filtered.loc[(cl_filtered['GIPM Z (OMNI mean)'] < 3) & (cl_filtered['GIPM Z (OMNI mean)'] > 1)]
cl_filtered_2 = cl_filtered.loc[(cl_filtered['GIPM Z (OMNI mean)'] < 5) & (cl_filtered['GIPM Z (OMNI mean)'] > 3)]
cl_filtered_3 = cl_filtered.loc[(cl_filtered['GIPM Z (OMNI mean)'] < 7) & (cl_filtered['GIPM Z (OMNI mean)'] > 9)]

##FIRST PLOT BATCH: Peak Frequency (Compressive & Transverse) and ellipticity split by MA/CA
MA_bounds = [5,10]
CA_bounds= [20,40]


Z_levels = {'Z+0': cl_filtered_0, 'Z+1': cl_filtered_1, 'Z+2': cl_filtered_2, 'Z+3': cl_filtered_3}

#filter by cone angle and mach no:

def ca_ma_filter_df_inp(df, ca_lims, ma_lims):
    ca_filt = df.loc[(df['cone angle (mean)'] >= ca_lims[0]) & (df['cone angle (mean)'] < ca_lims[1])]
    ca_ma_filt = ca_filt.loc[((ca_filt['M_A (mean)'] >= ma_lims[0]) & (ca_filt['M_A (mean)'] < ma_lims[1]))]
    return(ca_ma_filt)

filtered_Z_levels = {}

for df_key, df in Z_levels.items():
    filtered_Z_levels[df_key] = ca_ma_filter_df_inp(df, CA_bounds, MA_bounds)

print("dfs filtered")

#produce x and y bin edge lists, needed for later plots.

x_bin_edges = range(20)
y_bin_edges = range(-20, 20)

example_df = filtered_Z_levels["Z+0"]

_, xedg, yedg = np.histogram2d(
        example_df['GIPM X (OMNI mean)'].to_numpy(),
        example_df['GIPM Y (OMNI mean)'].to_numpy(),
        bins=[x_bin_edges, y_bin_edges]
    )

#heatmaps, with bins w/ under 30 obs removed.

histograms = {}
compressive_heatmap = {}
transverse_heatmap = {}

for group_name, level_df in filtered_Z_levels.items():
    histograms = {}
    compressive_heatmap = {}
    transverse_heatmap = {}
    histograms[group_name], compressive_heatmap[group_name], transverse_heatmap[group_name] = compute_hists2d_low_data(level_df)

print("ready to plot")

ThreeD_Binned_Plot(property_key, threeD_blocks, xedg, yedg)