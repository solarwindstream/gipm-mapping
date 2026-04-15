#Split location-only results by MA and cone angle, then go into those bins and pull out some spectra

import pandas as pd
import numpy as np
import datetime as dt
import glob

from XMA_finder import XMA_finder
from heatmap_modules import compute_freq_ellip_hists

import matplotlib
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as ticker

from heatmap_plotting_modules import draw_background, draw_hist, set_limits, mask_inside_magnetopause

##load Cluster CSVs

cl_file_list = []

path = "/Users/roseatkinson/Documents/Cluster_Integrated_CSVs/**"

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

#Filtering: remove times when OMNI source sc was too far from X-line, & current sheets, & filter down to -5<Zgipm<5 "plane"
cl_filtered = cl_power_all.loc[(cl_power_all['OMNI Dist from X line (mean)'] < 70) & (cl_power_all['Max IMF Deviation'] < 60)]
cl_filtered_lowZ = cl_filtered.loc[(cl_filtered['GIPM Z (OMNI mean)'] < 5) & (cl_filtered['GIPM Z (OMNI mean)'] > -5)]

##FIRST PLOT BATCH: Peak Frequency (Compressive & Transverse) and ellipticity split by MA/CA
MA_bounds = {"5_10": [5,10],"10_15": [10,15]}
CA_bounds_narrow = {"rad": [0,30], "lowspir": [30,45], "highspir": [45,60],"lowperp": [60,75], "highperp": [75,90]}

#filter by cone angle and mach no:

def ca_ma_filter(ca_lims, ma_lims):
    ca_filt = cl_filtered_lowZ.loc[(cl_filtered_lowZ['cone angle (mean)'] >= ca_lims[0]) & (cl_filtered_lowZ['cone angle (mean)'] < ca_lims[1])]
    ca_ma_filt = ca_filt.loc[((ca_filt['M_A (mean)'] >= ma_lims[0]) & (ca_filt['M_A (mean)'] < ma_lims[1]))]
    return(ca_ma_filt)

CA_MA_filtered_frames = {}

for ca_key, ca_bounds in CA_bounds.items:
    CA_MA_filtered_frames[ca_key] = {}
    for ma_key, ma_bounds in MA_bounds.items:
        CA_MA_filtered_frames[ca_key][ma_key] = ca_ma_filter(ca_bounds, ma_bounds)
    

#coverage histograms
x_bin_edges = range(20)
y_bin_edges = range(-20, 20)

#produce x and y bin edge lists, needed for later plots.
example_df = CA_MA_filtered_frames["rad"]["5_10"]

_, xedg, yedg = np.histogram2d(
        example_df['GIPM X (OMNI mean)'].to_numpy(),
        example_df['GIPM Y (OMNI mean)'].to_numpy(),
        bins=[x_bin_edges, y_bin_edges]
    )

#updated to include heatmaps, with bins w/ under 50 obs removed.

histograms = {}
peak_comp_freq_heatmap = {}
peak_trans_freq_heatmap = {}
ellipticity_heatmap = {}

for group_name, subsets in CA_MA_filtered_frames.items():
    histograms[group_name] = {}
    compressive_freq_heatmap[group_name] = {}
    transverse_freq_heatmap[group_name] = {}
    ellipticity_heatmap[group_name] = {}
    for subset_name, df in subsets.items():
        df = CA_MA_filtered_frames[group_name][subset_name]
        histograms[group_name][subset_name], peak_comp_freq_heatmap[group_name][subset_name], peak_trans_freq_heatmap[group_name][subset_name], ellipticity_heatmap_heatmap[group_name][subset_name] = compute_freq_ellip_hists(df)

#Make plots
#Save Plots!!