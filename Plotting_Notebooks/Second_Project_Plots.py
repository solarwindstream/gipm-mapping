import pandas as pd
import numpy as np
import datetime as dt
import glob
from heatmap_modules import compute_freq_ellip_hists, compute_normalised_freq_hists, compute_error_hists

from heatmap_plotting_modules import CA_MA_Binned_Plot, CA_Error_Plot

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

#Add a row with the Takahashi and Heilig predictions

#Takahashi prefactor
pi = np.pi
m_p = 1.67E-27
q_p = 1.60E-19
#1/2pi * q/m from gyrofreq
takahashi_pref = q_p/(4*pi*m_p)

cl_filtered_lowZ['Takahashi Frequency, Hz'] = 1E-9*takahashi_pref*(cl_filtered_lowZ['IMF B (mean)'])*(np.cos(np.deg2rad(cl_filtered_lowZ['cone angle (mean)'])))**2
cl_filtered_lowZ['Heilig Frequency, Hz'] = 1E-3*(0.78*cl_filtered_lowZ['M_A (mean)'] + 0.64)*cl_filtered_lowZ['IMF B (mean)']

cl_filtered_lowZ['Takahashi Transverse Error'] = (cl_filtered_lowZ['Takahashi Frequency, Hz'] - cl_filtered_lowZ['Peak Transverse Frequency'])/cl_filtered_lowZ['Peak Transverse Frequency']
cl_filtered_lowZ['Takahashi Compressive Error'] = (cl_filtered_lowZ['Takahashi Frequency, Hz'] - cl_filtered_lowZ['Peak Compressive Frequency'])/cl_filtered_lowZ['Peak Compressive Frequency']
cl_filtered_lowZ['Heilig Transverse Error'] = (cl_filtered_lowZ['Heilig Frequency, Hz'] - cl_filtered_lowZ['Peak Transverse Frequency'])/cl_filtered_lowZ['Peak Transverse Frequency']

print(cl_filtered_lowZ['Takahashi Transverse Error'].describe())
print(cl_filtered_lowZ['Takahashi Compressive Error'].describe())
print(cl_filtered_lowZ['Heilig Transverse Error'].describe())

##FIRST PLOT BATCH: Peak Frequency (Compressive & Transverse) and ellipticity split by MA/CA
MA_bounds = {"5_10": [5,10],"10_15": [10,15]}
CA_bounds_narrow = {"rad": [0,30], "lowspir": [30,45], "highspir": [45,60],"lowperp": [60,75], "highperp": [75,90]}

#filter by cone angle and mach no:

def ca_ma_filter(ca_lims, ma_lims):
    ca_filt = cl_filtered_lowZ.loc[(cl_filtered_lowZ['cone angle (mean)'] >= ca_lims[0]) & (cl_filtered_lowZ['cone angle (mean)'] < ca_lims[1])]
    ca_ma_filt = ca_filt.loc[((ca_filt['M_A (mean)'] >= ma_lims[0]) & (ca_filt['M_A (mean)'] < ma_lims[1]))]
    return(ca_ma_filt)

CA_MA_filtered_frames = {}

for ca_key, ca_bounds in CA_bounds_narrow.items():
    CA_MA_filtered_frames[ca_key] = {}
    for ma_key, ma_bounds in MA_bounds.items():
        CA_MA_filtered_frames[ca_key][ma_key] = ca_ma_filter(ca_bounds, ma_bounds)

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

# histograms = {}
# peak_comp_freq_heatmap = {}
# peak_trans_freq_heatmap = {}
# ellipticity_heatmap = {}

# for group_name, subsets in CA_MA_filtered_frames.items():
#     histograms[group_name] = {}
#     peak_comp_freq_heatmap[group_name] = {}
#     peak_trans_freq_heatmap[group_name] = {}
#     ellipticity_heatmap[group_name] = {}
#     for subset_name, df in subsets.items():
#         df = CA_MA_filtered_frames[group_name][subset_name]
#         histograms[group_name][subset_name], peak_comp_freq_heatmap[group_name][subset_name], peak_trans_freq_heatmap[group_name][subset_name], ellipticity_heatmap[group_name][subset_name] = compute_freq_ellip_hists(df)

#also generate Mach number difference and ratios for frequency and ellipticity

# trans_freq_ratios = {}
# trans_freq_diffs = {}

# for group_name, subsets in peak_trans_freq_heatmap.items():
#     trans_freq_ratios[group_name] = peak_trans_freq_heatmap[group_name]["10_15"]/peak_trans_freq_heatmap[group_name]["5_10"]
#     trans_freq_diffs[group_name] = peak_trans_freq_heatmap[group_name]["10_15"] - peak_trans_freq_heatmap[group_name]["5_10"]

# comp_freq_ratios = {}
# comp_freq_diffs = {}

# for group_name, subsets in peak_comp_freq_heatmap.items():
#     comp_freq_ratios[group_name] = peak_comp_freq_heatmap[group_name]["10_15"]/peak_comp_freq_heatmap[group_name]["5_10"]
#     comp_freq_diffs[group_name] = peak_comp_freq_heatmap[group_name]["10_15"] - peak_comp_freq_heatmap[group_name]["5_10"]

# ellipticity_ratios = {}
# ellipticity_diffs = {}

# for group_name, subsets in ellipticity_heatmap.items():
#     ellipticity_ratios[group_name] = ellipticity_heatmap[group_name]["10_15"]/ellipticity_heatmap[group_name]["5_10"]
#     ellipticity_diffs[group_name] = ellipticity_heatmap[group_name]["10_15"] - ellipticity_heatmap[group_name]["5_10"]

#now how to insert into plots? keep it basic for now

# trans_freq_blocks_ratio = [
#     [peak_trans_freq_heatmap["rad"]["10_15"], peak_trans_freq_heatmap["lowspir"]["10_15"], peak_trans_freq_heatmap["highspir"]["10_15"], peak_trans_freq_heatmap["lowperp"]["10_15"], peak_trans_freq_heatmap["highperp"]["10_15"]],
#     [peak_trans_freq_heatmap["rad"]["5_10"], peak_trans_freq_heatmap["lowspir"]["5_10"], peak_trans_freq_heatmap["highspir"]["5_10"], peak_trans_freq_heatmap["lowperp"]["5_10"], peak_trans_freq_heatmap["highperp"]["5_10"]],
#     [trans_freq_ratios["rad"], trans_freq_ratios["lowspir"], trans_freq_ratios["highspir"], trans_freq_ratios["lowperp"], trans_freq_ratios["highperp"]]
# ]

# trans_freq_blocks_diff = [
#     [peak_trans_freq_heatmap["rad"]["10_15"], peak_trans_freq_heatmap["lowspir"]["10_15"], peak_trans_freq_heatmap["highspir"]["10_15"], peak_trans_freq_heatmap["lowperp"]["10_15"], peak_trans_freq_heatmap["highperp"]["10_15"]],
#     [peak_trans_freq_heatmap["rad"]["5_10"], peak_trans_freq_heatmap["lowspir"]["5_10"], peak_trans_freq_heatmap["highspir"]["5_10"], peak_trans_freq_heatmap["lowperp"]["5_10"], peak_trans_freq_heatmap["highperp"]["5_10"]],
#     [trans_freq_diffs["rad"], trans_freq_diffs["lowspir"], trans_freq_diffs["highspir"], trans_freq_diffs["lowperp"], trans_freq_diffs["highperp"]]
# ]

# comp_freq_blocks_ratio = [
#     [peak_comp_freq_heatmap["rad"]["10_15"], peak_comp_freq_heatmap["lowspir"]["10_15"], peak_comp_freq_heatmap["highspir"]["10_15"], peak_comp_freq_heatmap["lowperp"]["10_15"], peak_comp_freq_heatmap["highperp"]["10_15"]],
#     [peak_comp_freq_heatmap["rad"]["5_10"], peak_comp_freq_heatmap["lowspir"]["5_10"], peak_comp_freq_heatmap["highspir"]["5_10"], peak_comp_freq_heatmap["lowperp"]["5_10"], peak_comp_freq_heatmap["highperp"]["5_10"]],
#     [comp_freq_ratios["rad"], comp_freq_ratios["lowspir"], comp_freq_ratios["highspir"], comp_freq_ratios["lowperp"], comp_freq_ratios["highperp"]]
# ]

# comp_freq_blocks_diff = [
#     [peak_comp_freq_heatmap["rad"]["10_15"], peak_comp_freq_heatmap["lowspir"]["10_15"], peak_comp_freq_heatmap["highspir"]["10_15"], peak_comp_freq_heatmap["lowperp"]["10_15"], peak_comp_freq_heatmap["highperp"]["10_15"]],
#     [peak_comp_freq_heatmap["rad"]["5_10"], peak_comp_freq_heatmap["lowspir"]["5_10"], peak_comp_freq_heatmap["highspir"]["5_10"], peak_comp_freq_heatmap["lowperp"]["5_10"], peak_comp_freq_heatmap["highperp"]["5_10"]],
#     [comp_freq_diffs["rad"], comp_freq_diffs["lowspir"], comp_freq_diffs["highspir"], comp_freq_diffs["lowperp"], comp_freq_diffs["highperp"]]
# ]

# ellipticity_blocks_ratio = [
#     [ellipticity_heatmap["rad"]["10_15"], ellipticity_heatmap["lowspir"]["10_15"], ellipticity_heatmap["highspir"]["10_15"], ellipticity_heatmap["lowperp"]["10_15"], ellipticity_heatmap["highperp"]["10_15"]],
#     [ellipticity_heatmap["rad"]["5_10"], ellipticity_heatmap["lowspir"]["5_10"], ellipticity_heatmap["highspir"]["5_10"], ellipticity_heatmap["lowperp"]["5_10"], ellipticity_heatmap["highperp"]["5_10"]],
#     [ellipticity_ratios["rad"], ellipticity_ratios["lowspir"], ellipticity_ratios["highspir"], ellipticity_ratios["lowperp"], ellipticity_ratios["highperp"]]
# ]

# ellipticity_blocks_diff = [
#     [ellipticity_heatmap["rad"]["10_15"], ellipticity_heatmap["lowspir"]["10_15"], ellipticity_heatmap["highspir"]["10_15"], ellipticity_heatmap["lowperp"]["10_15"], ellipticity_heatmap["highperp"]["10_15"]],
#     [ellipticity_heatmap["rad"]["5_10"], ellipticity_heatmap["lowspir"]["5_10"], ellipticity_heatmap["highspir"]["5_10"], ellipticity_heatmap["lowperp"]["5_10"], ellipticity_heatmap["highperp"]["5_10"]],
#     [ellipticity_diffs["rad"], ellipticity_diffs["lowspir"], ellipticity_diffs["highspir"], ellipticity_diffs["lowperp"], ellipticity_diffs["highperp"]]
# ]

#print("ready to plot")
#eventually do it all with dict but this will do for now.
#property keys: "Peak Transverse Frequency", "Peak Compressive Frequency", "Ellipticity","Normalised Compressive Frequency","Normalised Transverse Frequency"

# CA_MA_Binned_Plot("Peak Transverse Frequency", "Ratio", trans_freq_blocks_ratio, xedg, yedg)
# CA_MA_Binned_Plot("Peak Compressive Frequency", "Ratio", comp_freq_blocks_ratio, xedg, yedg)
# CA_MA_Binned_Plot("Ellipticity", "Ratio", ellipticity_blocks_ratio, xedg, yedg)

#now do histograms with NORMALISED frequencies!

#updated to include heatmaps, with bins w/ under 50 obs removed.

# histograms = {}
# norm_comp_freq_heatmap = {}
# norm_trans_freq_heatmap = {}

# for group_name, subsets in CA_MA_filtered_frames.items():
#     histograms[group_name] = {}
#     norm_comp_freq_heatmap[group_name] = {}
#     norm_trans_freq_heatmap[group_name] = {}
#     for subset_name, df in subsets.items():
#         df = CA_MA_filtered_frames[group_name][subset_name]
#         histograms[group_name][subset_name], norm_comp_freq_heatmap[group_name][subset_name], norm_trans_freq_heatmap[group_name][subset_name] = compute_normalised_freq_hists(df)

#also generate Mach number difference and ratios for frequency and ellipticity

# norm_trans_freq_ratios = {}
# norm_trans_freq_diffs = {}

# for group_name, subsets in norm_trans_freq_heatmap.items():
#     norm_trans_freq_ratios[group_name] = norm_trans_freq_heatmap[group_name]["10_15"]/norm_trans_freq_heatmap[group_name]["5_10"]
#     norm_trans_freq_diffs[group_name] = norm_trans_freq_heatmap[group_name]["10_15"] - norm_trans_freq_heatmap[group_name]["5_10"]

# norm_comp_freq_ratios = {}
# norm_comp_freq_diffs = {}

# for group_name, subsets in norm_comp_freq_heatmap.items():
#     norm_comp_freq_ratios[group_name] = norm_comp_freq_heatmap[group_name]["10_15"]/norm_comp_freq_heatmap[group_name]["5_10"]
#     norm_comp_freq_diffs[group_name] = norm_comp_freq_heatmap[group_name]["10_15"] - norm_comp_freq_heatmap[group_name]["5_10"]

# #now how to insert into plots? keep it basic for now

# norm_trans_freq_blocks_ratio = [
#     [norm_trans_freq_heatmap["rad"]["10_15"], norm_trans_freq_heatmap["lowspir"]["10_15"], norm_trans_freq_heatmap["highspir"]["10_15"], norm_trans_freq_heatmap["lowperp"]["10_15"], norm_trans_freq_heatmap["highperp"]["10_15"]],
#     [norm_trans_freq_heatmap["rad"]["5_10"], norm_trans_freq_heatmap["lowspir"]["5_10"], norm_trans_freq_heatmap["highspir"]["5_10"], norm_trans_freq_heatmap["lowperp"]["5_10"], norm_trans_freq_heatmap["highperp"]["5_10"]],
#     [norm_trans_freq_ratios["rad"], norm_trans_freq_ratios["lowspir"], norm_trans_freq_ratios["highspir"], norm_trans_freq_ratios["lowperp"], norm_trans_freq_ratios["highperp"]]
# ]

# norm_comp_freq_blocks_ratio = [
#     [norm_comp_freq_heatmap["rad"]["10_15"], norm_comp_freq_heatmap["lowspir"]["10_15"], norm_comp_freq_heatmap["highspir"]["10_15"], norm_comp_freq_heatmap["lowperp"]["10_15"], norm_comp_freq_heatmap["highperp"]["10_15"]],
#     [norm_comp_freq_heatmap["rad"]["5_10"], norm_comp_freq_heatmap["lowspir"]["5_10"], norm_comp_freq_heatmap["highspir"]["5_10"], norm_comp_freq_heatmap["lowperp"]["5_10"], norm_comp_freq_heatmap["highperp"]["5_10"]],
#     [comp_freq_ratios["rad"], comp_freq_ratios["lowspir"], comp_freq_ratios["highspir"], comp_freq_ratios["lowperp"], comp_freq_ratios["highperp"]]
# ]

# CA_MA_Binned_Plot("Normalised Transverse Frequency", "Ratio", norm_trans_freq_blocks_ratio, xedg, yedg)
# CA_MA_Binned_Plot("Normalised Compressive Frequency", "Ratio", norm_comp_freq_blocks_ratio, xedg, yedg)

#Takahashi and Heilig comparisons

#Split only by cone angle, not M_A

CA_bounds_narrow = {"rad": [0,30], "lowspir": [30,45], "highspir": [45,60],"lowperp": [60,75], "highperp": [75,90]}

#filter by cone angle and mach no:

def ca_filter(ca_lims):
    ca_filt = cl_filtered_lowZ.loc[(cl_filtered_lowZ['cone angle (mean)'] >= ca_lims[0]) & (cl_filtered_lowZ['cone angle (mean)'] < ca_lims[1])]
    return(ca_filt)

CA_filtered_frames = {}

for ca_key, ca_bounds in CA_bounds_narrow.items():
    CA_filtered_frames[ca_key] = ca_filter(ca_bounds)

print("CA-only filtered")

#updated to include heatmaps, with bins w/ under 50 obs removed.

histograms = []
takahashi_trans_error_heatmap = []
takahashi_comp_error_heatmap = []
heilig_trans_error_heatmap = []

for group_name, df in CA_filtered_frames.items():
    print(group_name)
    histogram, takahashi_trans_error_hm, takahashi_comp_error_hm, heilig_trans_error_hm = compute_error_hists(df)
    histograms.append(histogram)
    takahashi_trans_error_heatmap.append(takahashi_trans_error_hm)
    takahashi_comp_error_heatmap.append(takahashi_comp_error_hm)
    heilig_trans_error_heatmap.append(heilig_trans_error_hm)

print("histograms produced")

CA_Error_Plot('Takahashi Transverse Error', takahashi_trans_error_heatmap, xedg, yedg)
CA_Error_Plot('Takahashi Compressive Error', takahashi_comp_error_heatmap, xedg, yedg)
CA_Error_Plot('Heilig Transverse Error', heilig_trans_error_heatmap, xedg, yedg)