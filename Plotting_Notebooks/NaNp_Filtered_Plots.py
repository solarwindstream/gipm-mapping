#Split location-only results by MA and cone angle, then go into those bins and pull out some spectra

import pandas as pd
import numpy as np
import datetime as dt
import glob
from heatmap_modules import compute_hists2d, compute_freq_ellip_hists, compute_normalised_freq_hists, compute_freq_hists_low_data

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

from heatmap_plotting_modules import model_calcs, cone_angle_line, draw_background, draw_heatmap, set_limits, mask_inside_magnetopause, NaNp_heatmap_plot, NaNp_heatmap_plot_B_filt

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

#add a row normalised by OMNI mag field average... I should be doing this at the source really.
cl_filtered_lowZ['Normalised Transverse Frequency'] = (cl_filtered_lowZ['Peak Transverse Frequency']/cl_filtered_lowZ['IMF B (mean)'])
cl_filtered_lowZ['Normalised Compressive Frequency'] = (cl_filtered_lowZ['Peak Compressive Frequency']/cl_filtered_lowZ['IMF B (mean)'])

#temporary B filter:

cl_filtered_midB = cl_filtered_lowZ.loc[(cl_filtered_lowZ['IMF B (mean)'] > 3) & (cl_filtered_lowZ['IMF B (mean)'] < 5)]

#produce histogram of NaNp values for whole (filtered) dataset


# Plotting a basic histogram
plt.hist(cl_filtered_lowZ['SW Na/Np (mean)'], bins=[0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1], color='skyblue', edgecolor='black')

# Adding labels and title
plt.xlabel('NaNp')
plt.ylabel('Observations')
plt.xlim(0, 0.1)

path = "/Users/roseatkinson/Documents/New_Figs/Na_Np_Mean_Hist.png"
plt.savefig(path)


#now filter for 30-52.5deg, 52.5-75deg, and Na/Np <0.02, >0.03.

##FIRST PLOT BATCH: Peak Frequency (Compressive & Transverse)
NaNp_bounds = {"NaNp < 0.03": [0,0.03],"NaNp > 0.045": [0.045,1]}
CA_bounds_wide = {"30-52.5°": [30,52.5], "52.5-75°": [52.5,75]}

#filter by cone angle and mach no:

def ca_nanp_filter(ca_lims, nanp_lims):
    ca_filt = cl_filtered_lowZ.loc[(cl_filtered_lowZ['cone angle (mean)'] >= ca_lims[0]) & (cl_filtered_lowZ['cone angle (mean)'] < ca_lims[1])]
    ca_ma_filt = ca_filt.loc[((ca_filt['SW Na/Np (mean)'] >= nanp_lims[0]) & (ca_filt['SW Na/Np (mean)'] < nanp_lims[1]))]
    return(ca_ma_filt)

CA_NaNp_filtered_frames = {}

for ca_key, ca_bounds in CA_bounds_wide.items():
    CA_NaNp_filtered_frames[ca_key] = {}
    for na_key, na_bounds in NaNp_bounds.items():
        CA_NaNp_filtered_frames[ca_key][na_key] = ca_nanp_filter(ca_bounds, na_bounds)
        print(ca_key, na_key, 'shape:', CA_NaNp_filtered_frames[ca_key][na_key].shape)

print("dfs filtered")



#produce x and y bin edge lists, needed for later plots.

x_bin_edges = range(0, 20)
y_bin_edges = range(-20, 20)

example_df = CA_NaNp_filtered_frames["30-52.5°"]["NaNp < 0.03"]

_, xedg, yedg = np.histogram2d(
        example_df['GIPM X (OMNI mean)'].to_numpy(),
        example_df['GIPM Y (OMNI mean)'].to_numpy(),
        bins=[x_bin_edges, y_bin_edges]
    )

#frequency

histograms = {}
peak_comp_freq_heatmap = {}
peak_trans_freq_heatmap = {}

for group_name, subsets in CA_NaNp_filtered_frames.items():
    histograms[group_name] = {}
    peak_comp_freq_heatmap[group_name] = {}
    peak_trans_freq_heatmap[group_name] = {}
    for subset_name, df in subsets.items():
        df = CA_NaNp_filtered_frames[group_name][subset_name]
        histograms[group_name][subset_name], peak_comp_freq_heatmap[group_name][subset_name], peak_trans_freq_heatmap[group_name][subset_name] = compute_freq_hists_low_data(df)

#now need to produce ratios

#also generate Mach number difference and ratios for frequency and ellipticity

trans_freq_ratios = {}
trans_freq_diffs = {}

for group_name, subsets in peak_trans_freq_heatmap.items():
    trans_freq_ratios[group_name] = peak_trans_freq_heatmap[group_name]["NaNp > 0.045"]/peak_trans_freq_heatmap[group_name]["NaNp < 0.03"]
    trans_freq_diffs[group_name] = peak_trans_freq_heatmap[group_name]["NaNp > 0.045"] - peak_trans_freq_heatmap[group_name]["NaNp < 0.03"]

comp_freq_ratios = {}
comp_freq_diffs = {}

for group_name, subsets in peak_comp_freq_heatmap.items():
    comp_freq_ratios[group_name] = peak_comp_freq_heatmap[group_name]["NaNp > 0.045"]/peak_comp_freq_heatmap[group_name]["NaNp < 0.03"]
    comp_freq_diffs[group_name] = peak_comp_freq_heatmap[group_name]["NaNp > 0.045"] - peak_comp_freq_heatmap[group_name]["NaNp < 0.03"]
    
#now how to insert into plots? keep it basic for now

trans_freq_blocks_ratio = [
    [peak_trans_freq_heatmap["30-52.5°"]["NaNp > 0.045"], peak_trans_freq_heatmap["52.5-75°"]["NaNp > 0.045"]] ,
    [peak_trans_freq_heatmap["30-52.5°"]["NaNp < 0.03"], peak_trans_freq_heatmap["52.5-75°"]["NaNp < 0.03"]],
    [trans_freq_ratios["30-52.5°"], trans_freq_ratios["52.5-75°"]]
]

comp_freq_blocks_ratio = [
    [peak_comp_freq_heatmap["30-52.5°"]["NaNp > 0.045"], peak_comp_freq_heatmap["52.5-75°"]["NaNp > 0.045"]],
    [peak_comp_freq_heatmap["30-52.5°"]["NaNp < 0.03"], peak_comp_freq_heatmap["52.5-75°"]["NaNp < 0.03"]],
    [comp_freq_ratios["30-52.5°"], comp_freq_ratios["52.5-75°"]]
]

#now make new blocks just for high and low NaNp, and ratios

wide_angle_blocks = [
    [peak_trans_freq_heatmap["30-52.5°"]["NaNp > 0.045"], peak_trans_freq_heatmap["52.5-75°"]["NaNp > 0.045"], peak_comp_freq_heatmap["30-52.5°"]["NaNp > 0.045"], peak_comp_freq_heatmap["52.5-75°"]["NaNp > 0.045"]],
    [peak_trans_freq_heatmap["30-52.5°"]["NaNp < 0.03"], peak_trans_freq_heatmap["52.5-75°"]["NaNp < 0.03"], peak_comp_freq_heatmap["30-52.5°"]["NaNp < 0.03"], peak_comp_freq_heatmap["52.5-75°"]["NaNp < 0.03"]],
    [trans_freq_ratios['30-52.5°'], trans_freq_ratios["52.5-75°"], comp_freq_ratios['30-52.5°'], comp_freq_ratios["52.5-75°"]]
]


###NORMALISED PLOTS NOW!

histograms = {}
norm_comp_freq_heatmap = {}
norm_trans_freq_heatmap = {}

for group_name, subsets in CA_NaNp_filtered_frames.items():
    histograms[group_name] = {}
    norm_comp_freq_heatmap[group_name] = {}
    norm_trans_freq_heatmap[group_name] = {}
    for subset_name, df in subsets.items():
        df = CA_NaNp_filtered_frames[group_name][subset_name]
        histograms[group_name][subset_name], norm_comp_freq_heatmap[group_name][subset_name], norm_trans_freq_heatmap[group_name][subset_name] = compute_normalised_freq_hists(df)

#now need to produce ratios

norm_trans_freq_ratios = {}
norm_trans_freq_diffs = {}

for group_name, subsets in norm_trans_freq_heatmap.items():
    norm_trans_freq_ratios[group_name] = norm_trans_freq_heatmap[group_name]["NaNp > 0.045"]/norm_trans_freq_heatmap[group_name]["NaNp < 0.03"]
    norm_trans_freq_diffs[group_name] = norm_trans_freq_heatmap[group_name]["NaNp > 0.045"] - norm_trans_freq_heatmap[group_name]["NaNp < 0.03"]

norm_comp_freq_ratios = {}
norm_comp_freq_diffs = {}

for group_name, subsets in norm_comp_freq_heatmap.items():
    norm_comp_freq_ratios[group_name] = norm_comp_freq_heatmap[group_name]["NaNp > 0.045"]/norm_comp_freq_heatmap[group_name]["NaNp < 0.03"]
    norm_comp_freq_diffs[group_name] = norm_comp_freq_heatmap[group_name]["NaNp > 0.045"] - norm_comp_freq_heatmap[group_name]["NaNp < 0.03"]
    
#now how to insert into plots? keep it basic for now

norm_trans_freq_blocks_ratio = [
    [norm_trans_freq_heatmap["30-52.5°"]["NaNp > 0.045"], norm_trans_freq_heatmap["52.5-75°"]["NaNp > 0.045"]] ,
    [norm_trans_freq_heatmap["30-52.5°"]["NaNp < 0.03"], norm_trans_freq_heatmap["52.5-75°"]["NaNp < 0.03"]],
    [norm_trans_freq_ratios["30-52.5°"], norm_trans_freq_ratios["52.5-75°"]]
]

norm_comp_freq_blocks_ratio = [
    [norm_comp_freq_heatmap["30-52.5°"]["NaNp > 0.045"], norm_comp_freq_heatmap["52.5-75°"]["NaNp > 0.045"]],
    [norm_comp_freq_heatmap["30-52.5°"]["NaNp < 0.03"], norm_comp_freq_heatmap["52.5-75°"]["NaNp < 0.03"]],
    [norm_comp_freq_ratios["30-52.5°"], norm_comp_freq_ratios["52.5-75°"]]
]

#now make new blocks just for high and low NaNp, and ratios

norm_wide_angle_blocks = [
    [norm_trans_freq_heatmap["30-52.5°"]["NaNp > 0.045"], norm_trans_freq_heatmap["52.5-75°"]["NaNp > 0.045"], norm_comp_freq_heatmap["30-52.5°"]["NaNp > 0.045"], norm_comp_freq_heatmap["52.5-75°"]["NaNp > 0.045"]],
    [norm_trans_freq_heatmap["30-52.5°"]["NaNp < 0.03"], norm_trans_freq_heatmap["52.5-75°"]["NaNp < 0.03"], norm_comp_freq_heatmap["30-52.5°"]["NaNp < 0.03"], norm_comp_freq_heatmap["52.5-75°"]["NaNp < 0.03"]],
    [norm_trans_freq_ratios['30-52.5°'], norm_trans_freq_ratios["52.5-75°"], norm_comp_freq_ratios['30-52.5°'], norm_comp_freq_ratios["52.5-75°"]]
]

#print(peak_trans_freq_heatmap["30-52.5°"]["NaNp > 0.03"])

#"Peak Frequency" and "Normalised Frequency"

#NaNp_heatmap_plot("Peak Frequency", wide_angle_blocks, xedg, yedg, ['0.035', '0.045'])
#NaNp_heatmap_plot("Normalised Frequency", norm_wide_angle_blocks, xedg, yedg, ['0.035', '0.045'])
NaNp_heatmap_plot_B_filt("Peak Frequency", wide_angle_blocks, xedg, yedg, ['0.03', '0.045'])

