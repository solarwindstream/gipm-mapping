#Split location-only results by MA and cone angle, then go into those bins and pull out some spectra

import pandas as pd
import numpy as np
import datetime as dt
import glob
from heatmap_modules import compute_hists2d, compute_freq_ellip_hists

import matplotlib
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Rectangle
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import TransformedBbox
from matplotlib.collections import LineCollection
import matplotlib.collections as mcoll
import matplotlib.path as mpath

from heatmap_plotting_modules import model_calcs, cone_angle_line, draw_background, draw_heatmap, set_limits, mask_inside_magnetopause

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

#Filtering: remove times when OMNI source sc was too far from X-line, & current sheets, & filter down to -5<Zgipm<5 "plane"
cl_filtered = cl_power_all.loc[(cl_power_all['OMNI Dist from X line (mean)'] < 70) & (cl_power_all['Max IMF Deviation'] < 60)]
cl_filtered_lowZ = cl_filtered.loc[(cl_filtered['GIPM Z (OMNI mean)'] < 5) & (cl_filtered['GIPM Z (OMNI mean)'] > -5)]

#produce histogram of NaNp values for whole (filtered) dataset


# Plotting a basic histogram
plt.hist(cl_filtered_lowZ['SW Na/Np (mean)'], bins=[0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2], color='skyblue', edgecolor='black')

# Adding labels and title
plt.xlabel('NaNp')
plt.ylabel('Observations')

path = "/Users/roseatkinson/Documents/New_Figs/Na_Np_Mean_Hist.png"
plt.savefig(path)


#now filter for 30-52.5deg, 52.5-75deg, and Na/Np <0.01, >0.04.

##FIRST PLOT BATCH: Peak Frequency (Compressive & Transverse) and ellipticity split by MA/CA
NaNp_bounds = {"NaNp < 0.01": [0,0.01],"NaNp > 0.04": [0.04,1]}
CA_bounds_wide = {"30-52.5°": [30,52.5], "52.5-75°": [52.5,75]}

#filter by cone angle and mach no:

def ca_nanp_filter(ca_lims, nanp_lims):
    ca_filt = cl_filtered_lowZ.loc[(cl_filtered_lowZ['cone angle (mean)'] >= ca_lims[0]) & (cl_filtered_lowZ['cone angle (mean)'] < ca_lims[1])]
    ca_ma_filt = ca_filt.loc[((ca_filt['M_A (mean)'] >= nanp_lims[0]) & (ca_filt['M_A (mean)'] < nanp_lims[1]))]
    return(ca_ma_filt)

CA_NaNp_filtered_frames = {}

for ca_key, ca_bounds in CA_bounds_wide.items():
    CA_NaNp_filtered_frames[ca_key] = {}
    for na_key, na_bounds in NaNp_bounds.items():
        CA_NaNp_filtered_frames[ca_key][na_key] = ca_nanp_filter(ca_bounds, na_bounds)

print("dfs filtered")

#produce x and y bin edge lists, needed for later plots.

x_bin_edges = range(20)
y_bin_edges = range(-20, 20)

example_df = CA_MA_filtered_frames["rad"]["5_10"]

_, xedg, yedg = np.histogram2d(
        example_df['GIPM X (OMNI mean)'].to_numpy(),
        example_df['GIPM Y (OMNI mean)'].to_numpy(),
        bins=[x_bin_edges, y_bin_edges]
    )

#updated to include heatmaps, with bins w/ under 50 obs removed.

#power first

histograms = {}
comp_power_heatmap = {}
trans_power_heatmap = {}
compressibility_heatmap = {}

for group_name, subsets in CA_NaNp_filtered_frames.items():
    histograms[group_name] = {}
    comp_power_heatmap[group_name] = {}
    trans_power_heatmap[group_name] = {}
    ellipticity_heatmap[group_name] = {}
    for subset_name, df in subsets.items():
        df = CA_MA_filtered_frames[group_name][subset_name]
        histograms[group_name][subset_name], comp_power_heatmap[group_name][subset_name], trans_power_heatmap[group_name][subset_name], compressibility_heatmap[group_name][subset_name] = compute_hists2d(df)


#then frequency

histograms = {}
peak_comp_freq_heatmap = {}
peak_trans_freq_heatmap = {}
ellipticity_heatmap = {}

for group_name, subsets in CA_NaNp_filtered_frames.items():
    histograms[group_name] = {}
    peak_comp_freq_heatmap[group_name] = {}
    peak_trans_freq_heatmap[group_name] = {}
    ellipticity_heatmap[group_name] = {}
    for subset_name, df in subsets.items():
        df = CA_MA_filtered_frames[group_name][subset_name]
        histograms[group_name][subset_name], peak_comp_freq_heatmap[group_name][subset_name], peak_trans_freq_heatmap[group_name][subset_name], ellipticity_heatmap[group_name][subset_name] = compute_freq_ellip_hists(df)
