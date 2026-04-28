#Split location-only results by MA and cone angle, then go into those bins and pull out some spectra

import pandas as pd
import numpy as np
import datetime as dt
import glob
from heatmap_modules import compute_freq_ellip_hists

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

histograms = {}
peak_comp_freq_heatmap = {}
peak_trans_freq_heatmap = {}
ellipticity_heatmap = {}

for group_name, subsets in CA_MA_filtered_frames.items():
    histograms[group_name] = {}
    peak_comp_freq_heatmap[group_name] = {}
    peak_trans_freq_heatmap[group_name] = {}
    ellipticity_heatmap[group_name] = {}
    for subset_name, df in subsets.items():
        df = CA_MA_filtered_frames[group_name][subset_name]
        histograms[group_name][subset_name], peak_comp_freq_heatmap[group_name][subset_name], peak_trans_freq_heatmap[group_name][subset_name], ellipticity_heatmap[group_name][subset_name] = compute_freq_ellip_hists(df)

#also generate Mach number difference and ratios for frequency and ellipticity

trans_freq_ratios = {}
trans_freq_diffs = {}

for group_name, subsets in peak_trans_freq_heatmap.items():
    trans_freq_ratios[group_name] = peak_trans_freq_heatmap[group_name]["10_15"]/peak_trans_freq_heatmap[group_name]["5_10"]
    trans_freq_diffs[group_name] = peak_trans_freq_heatmap[group_name]["10_15"] - peak_trans_freq_heatmap[group_name]["5_10"]

comp_freq_ratios = {}
comp_freq_diffs = {}

for group_name, subsets in peak_comp_freq_heatmap.items():
    comp_freq_ratios[group_name] = peak_comp_freq_heatmap[group_name]["10_15"]/peak_comp_freq_heatmap[group_name]["5_10"]
    comp_freq_diffs[group_name] = peak_comp_freq_heatmap[group_name]["10_15"] - peak_comp_freq_heatmap[group_name]["5_10"]

ellipticity_ratios = {}
ellipticity_diffs = {}

for group_name, subsets in ellipticity_heatmap.items():
    ellipticity_ratios[group_name] = ellipticity_heatmap[group_name]["10_15"]/ellipticity_heatmap[group_name]["5_10"]
    ellipticity_diffs[group_name] = ellipticity_heatmap[group_name]["10_15"] - ellipticity_heatmap[group_name]["5_10"]

#now how to insert into plots? keep it basic for now

trans_freq_blocks_ratio = [
    [peak_trans_freq_heatmap["rad"]["10_15"], peak_trans_freq_heatmap["lowspir"]["10_15"], peak_trans_freq_heatmap["highspir"]["10_15"], peak_trans_freq_heatmap["lowperp"]["10_15"], peak_trans_freq_heatmap["highperp"]["10_15"]],
    [peak_trans_freq_heatmap["rad"]["5_10"], peak_trans_freq_heatmap["lowspir"]["5_10"], peak_trans_freq_heatmap["highspir"]["5_10"], peak_trans_freq_heatmap["lowperp"]["5_10"], peak_trans_freq_heatmap["highperp"]["5_10"]],
    [trans_freq_ratios["rad"], trans_freq_ratios["lowspir"], trans_freq_ratios["highspir"], trans_freq_ratios["lowperp"], trans_freq_ratios["highperp"]]
]

trans_freq_blocks_diff = [
    [peak_trans_freq_heatmap["rad"]["10_15"], peak_trans_freq_heatmap["lowspir"]["10_15"], peak_trans_freq_heatmap["highspir"]["10_15"], peak_trans_freq_heatmap["lowperp"]["10_15"], peak_trans_freq_heatmap["highperp"]["10_15"]],
    [peak_trans_freq_heatmap["rad"]["5_10"], peak_trans_freq_heatmap["lowspir"]["5_10"], peak_trans_freq_heatmap["highspir"]["5_10"], peak_trans_freq_heatmap["lowperp"]["5_10"], peak_trans_freq_heatmap["highperp"]["5_10"]],
    [trans_freq_diffs["rad"], trans_freq_diffs["lowspir"], trans_freq_diffs["highspir"], trans_freq_diffs["lowperp"], trans_freq_diffs["highperp"]]
]

comp_freq_blocks_ratio = [
    [peak_comp_freq_heatmap["rad"]["10_15"], peak_comp_freq_heatmap["lowspir"]["10_15"], peak_comp_freq_heatmap["highspir"]["10_15"], peak_comp_freq_heatmap["lowperp"]["10_15"], peak_comp_freq_heatmap["highperp"]["10_15"]],
    [peak_comp_freq_heatmap["rad"]["5_10"], peak_comp_freq_heatmap["lowspir"]["5_10"], peak_comp_freq_heatmap["highspir"]["5_10"], peak_comp_freq_heatmap["lowperp"]["5_10"], peak_comp_freq_heatmap["highperp"]["5_10"]],
    [comp_freq_ratios["rad"], comp_freq_ratios["lowspir"], comp_freq_ratios["highspir"], comp_freq_ratios["lowperp"], comp_freq_ratios["highperp"]]
]

comp_freq_blocks_diff = [
    [peak_comp_freq_heatmap["rad"]["10_15"], peak_comp_freq_heatmap["lowspir"]["10_15"], peak_comp_freq_heatmap["highspir"]["10_15"], peak_comp_freq_heatmap["lowperp"]["10_15"], peak_comp_freq_heatmap["highperp"]["10_15"]],
    [peak_comp_freq_heatmap["rad"]["5_10"], peak_comp_freq_heatmap["lowspir"]["5_10"], peak_comp_freq_heatmap["highspir"]["5_10"], peak_comp_freq_heatmap["lowperp"]["5_10"], peak_comp_freq_heatmap["highperp"]["5_10"]],
    [comp_freq_diffs["rad"], comp_freq_diffs["lowspir"], comp_freq_diffs["highspir"], comp_freq_diffs["lowperp"], comp_freq_diffs["highperp"]]
]

ellipticity_blocks_ratio = [
    [ellipticity_heatmap["rad"]["10_15"], ellipticity_heatmap["lowspir"]["10_15"], ellipticity_heatmap["highspir"]["10_15"], ellipticity_heatmap["lowperp"]["10_15"], ellipticity_heatmap["highperp"]["10_15"]],
    [ellipticity_heatmap["rad"]["5_10"], ellipticity_heatmap["lowspir"]["5_10"], ellipticity_heatmap["highspir"]["5_10"], ellipticity_heatmap["lowperp"]["5_10"], ellipticity_heatmap["highperp"]["5_10"]],
    [ellipticity_ratios["rad"], ellipticity_ratios["lowspir"], ellipticity_ratios["highspir"], ellipticity_ratios["lowperp"], ellipticity_ratios["highperp"]]
]

ellipticity_blocks_diff = [
    [ellipticity_heatmap["rad"]["10_15"], ellipticity_heatmap["lowspir"]["10_15"], ellipticity_heatmap["highspir"]["10_15"], ellipticity_heatmap["lowperp"]["10_15"], ellipticity_heatmap["highperp"]["10_15"]],
    [ellipticity_heatmap["rad"]["5_10"], ellipticity_heatmap["lowspir"]["5_10"], ellipticity_heatmap["highspir"]["5_10"], ellipticity_heatmap["lowperp"]["5_10"], ellipticity_heatmap["highperp"]["5_10"]],
    [ellipticity_diffs["rad"], ellipticity_diffs["lowspir"], ellipticity_diffs["highspir"], ellipticity_diffs["lowperp"], ellipticity_diffs["highperp"]]
]

print("ready to plot")
#eventually do it all with dict but this will do for now.
#property keys: "Peak Transverse Frequency", "Peak Compressive Frequency", "Ellipticity"

def CA_MA_Binned_Plot(property_key,comparison_key, ma_ca_blocks):
    """Cone Angle and Mach Number 3x5 grid """
    #find constants for Merka & Shue
    X_shue, R_shue, Xgipm, Ygipm, Zgipm, f, fitting_coeffs = model_calcs()

    #main cbar titles:

    cbar_titles = {"Peak Transverse Frequency": "Frequency, Hz", "Peak Compressive Frequency": "Frequency, Hz", "Ellipticity": "Ellipticity"}

    fig = plt.figure(figsize=(9, 8), dpi=300, constrained_layout=True)
    gs = fig.add_gridspec(
        nrows=3, ncols=6,      # 1 column for patch labels
        width_ratios=[0.35, 1, 1, 1, 1, 1],  # label column thinner
        wspace=0.05, hspace=0.1
    )

    fig.suptitle(property_key, fontsize=18)
    plt.rcParams['axes.labelsize'] = 14

    # Row labels (top row → bottom row)
    row_labels = [
        r'$10 \leq M_A < 15$',
        r'$5 \leq M_A < 10$',
        comparison_key
    ]
    angle_titles = ["0–30°", "30–45°", "45–60°", "60–75°", "75–90°"]
    extent = [xedg[0], xedg[-1], yedg[0], yedg[-1]]
    # -------------------------------
    # MAKE AXES FOR THE 3×5 PANELS
    # -------------------------------
    axs = []
    for r in range(3):
        plot_row_axes = []
        gs_row = r      
        for c in range(5):
            ax = fig.add_subplot(gs[gs_row, c + 1])
            plot_row_axes.append(ax)
        axs.append(plot_row_axes)

    # -------------------------------
    # Patch Labels (Rounded Boxes)
    # -------------------------------

    for r in range(3):
        ax_patch = fig.add_subplot(gs[r, 0])
        ax_patch.set_axis_off()

        # -- Draw text first so we can query its bounding box --
        txt = ax_patch.text(
            0.5, 0.5,                     # centered in the Axes
            row_labels[r],
            ha="center",
            va="center",
            fontsize=12,
            transform=ax_patch.transAxes,
            rotation='vertical'
        )

        fig.canvas.draw()  # required to obtain correct text bounding box

        # -- Convert text bounding box from display to Axes coordinates --
        renderer = fig.canvas.get_renderer()
        bbox = txt.get_window_extent(renderer=renderer)
        bbox_axes = TransformedBbox(
            bbox, ax_patch.transAxes.inverted()
        )

        # Add some padding around the text
        pad_x = 0.04   # fractional padding in axes coordinates
        pad_y = 0.01

        x0 = bbox_axes.x0 - pad_x
        y0 = bbox_axes.y0 - pad_y
        width = bbox_axes.width + 2 * pad_x
        height = bbox_axes.height + 2 * pad_y

        # -- Rounded box placed behind the text --
        box = FancyBboxPatch(
            (x0, y0),
            width,
            height,
            boxstyle="round,pad=0.2,rounding_size=0.06",
            fc="lightgrey",
            ec="dimgrey",
            linewidth=1,
            mutation_aspect=1,
            transform=ax_patch.transAxes,
            zorder=0.5,
        )
        ax_patch.add_patch(box)

        # Move text above box
        txt.set_zorder(1)


    # -------------------------------
    # COLORMAP
    # -------------------------------

    cmap_dict = {"Peak Transverse Frequency": [0.001, 0.1], "Peak Compressive Frequency": [0.001, 0.1], "Ellipticity": [1,4]}
    powercmp = 'magma'
    #make norm, not lognorm 
    if property_key=='Ellipticity':
        power_norm = colors.Normalize(vmin=cmap_dict[property_key][0], vmax=cmap_dict[property_key][1])
    else:
        power_norm = colors.LogNorm(vmin=cmap_dict[property_key][0], vmax=cmap_dict[property_key][1])
    
    comparison_cmp = 'RdBu_r'
    if comparison_key == "Ratio":
        comparison_norm = colors.LogNorm(vmin=0.1, vmax=10)
    else:
        comparison_norm = colors.LogNorm(vmin=0.1, vmax=10)

    # -------------------------------
    # PLOT ALL PANELS
    # -------------------------------

    for col in range(5):    # angle class
        
        title = angle_titles[col]

        for row in range(3):                     # mach no. class
            ax = axs[row][col]

            # Draw contour, magnetopause
            draw_background(ax, Xgipm[:, :, 0], Ygipm[:, :, 0], f[:, :, 0],
                            X_shue, R_shue)

            # Histogram for this cell
            hist = ma_ca_blocks[row][col]
            
            angle_line = cone_angle_line(fitting_coeffs, title)
            
            if row < 2:
                draw_heatmap(ax, hist, extent, powercmp, power_norm, angle_line)

            if row == 2:
                draw_heatmap(ax, hist, extent, comparison_cmp, comparison_norm, angle_line)
                
            mask_inside_magnetopause(ax, X_shue, R_shue)
            
            # redraw magnetopause boundary so it stays crisp
            ax.plot(X_shue, R_shue, 'k', linewidth=1)
            set_limits(ax)

            # Labels
            if col == 0:
                ax.set_ylabel("$Y_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")
            if row == 0:
                ax.set_title(rf'$\alpha$ = {title}', fontsize=12)
            if row == 2:
                ax.set_xlabel("$X_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")

    # -------------------------------
    # COLORBARS (TWO SEPARATE ON RIGHT)
    # -------------------------------

    from matplotlib.cm import ScalarMappable

    # --- Scalar mappables (independent of any single subplot image)
    sm_power = ScalarMappable(norm=power_norm, cmap=powercmp)
    sm_power.set_array([])

    sm_ratio = ScalarMappable(norm=comparison_norm, cmap=comparison_cmp)
    sm_ratio.set_array([])

    # --- Top two rows colourbar
    top_axes = axs[0] + axs[1]   # flatten row 0 and 1
    cbar1 = fig.colorbar(
        sm_power,
        ax=top_axes,
        location='right',
        pad=0.02,
        extend='both'
    )
    cbar1.set_label(cbar_titles[property_key])

    # --- Bottom row colourbar (ratio)
    bottom_axes = axs[2]
    cbar2 = fig.colorbar(
        sm_ratio,
        ax=bottom_axes,
        location='right',
        aspect=10,
        pad=0.02,
        extend='both'
    )
    cbar2.set_label(comparison_key)

    cbar1.ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext())
    cbar2.ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext())

    path = "/Users/roseatkinson/Documents/New_Figs/CA_MA_" + property_key + comparison_key + ".png"
    plt.savefig(path)
    #Save Plots!!

CA_MA_Binned_Plot("Peak Transverse Frequency", "Ratio", trans_freq_blocks_ratio)
CA_MA_Binned_Plot("Peak Compressive Frequency", "Ratio", comp_freq_blocks_ratio)
CA_MA_Binned_Plot("Ellipticity", "Ratio", ellipticity_blocks_ratio)