import pandas as pd
import numpy as np
import datetime as dt

from merka05_surface_eq_array_GIPM import merka05_surface_eq_array_GIPM

import matplotlib
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import TransformedBbox
from matplotlib.collections import LineCollection
import matplotlib.collections as mcoll
import matplotlib.path as mpath
from cmcrameri import cm as cm_cram

def model_calcs():
    """Calculate parameters for Shue and Merka models required for plot"""

    #Shue magnetopause model, D = 1.76 nPa

    m_1 = 10.22
    m_2 = 1.29
    m_3 = 0.184
    m_4 = 8.14
    m_5 = 6.6
    m_6 = 0.58
    m_7 = -0.007
    m_8 = 0.024

    #use B_z = 0
    B_z = 0
    D_p = 1.76

    alpha = (m_6 + m_7*B_z)*(1 +m_8*(np.log(D_p)))

    tanh_angle = m_3*(B_z+m_4)
    r_0 = (m_1 + m_2*np.tanh(tanh_angle))*(D_p**(-1/m_5))

    pi = np.pi

    theta = np.arange(-pi/2, pi/2, 0.01)

    r_mod = (2/(1+np.cos(theta)))**alpha
    r = r_0*r_mod

    X_shue = r*(np.cos(theta))
    R_shue = r*(np.sin(theta))

    #make a coverage plot

    x = np.linspace(0, 20, 100) #x coordinates (Re)
    y = np.linspace(-30, 30, 100) #y coordinates (Re)
    z = 0 #z coordinates in Re

    [Xgipm,Ygipm,Zgipm] = np.meshgrid(x,y,z,indexing="ij")
    XMA_all = 10
    fitting_coeffs = merka05_surface_eq_array_GIPM(XMA_all)

    Xn = Xgipm
    Yn = Ygipm
    Zn = Zgipm
    f = fitting_coeffs[0]*Xn**2 + fitting_coeffs[1]*Yn**2 + fitting_coeffs[2]*Zn**2+ 2*fitting_coeffs[3]*Xn*Yn + 2*fitting_coeffs[4]*Yn*Zn + 2*fitting_coeffs[5]*Xn*Zn + 2*fitting_coeffs[6]*Xn+2*fitting_coeffs[7]*Yn + 2*fitting_coeffs[8]*Zn + fitting_coeffs[9]
    
    return(X_shue, R_shue, Xn, Yn, Zn, f, fitting_coeffs)


def draw_background(ax, xg, yg, f, x_shue, r_shue):
    """Draw bow shock, magnetopause, and y=0 line"""
    ax.contour(xg, yg, f, levels=[0], colors="black", linewidths=1)
    ax.plot(x_shue, r_shue, 'k', linewidth=1)
    ax.hlines(0, 0, 25, color='k', linewidth=1)

def cone_angle_line(fitting_coeffs, angle_key):
    # Line slopes for different angle classes
    line_slopes = {
        "0–30°": np.tan(np.deg2rad(15)),
        "30–45°": np.tan(np.deg2rad(37.5)),
        "45–60°": np.tan(np.deg2rad(52.5)),
        "60–75°": np.tan(np.deg2rad(67.5)),
        "75–90°": np.tan(np.deg2rad(82.5)),
        "30-52.5°": np.tan(np.deg2rad(41.25)),
        "52.5-75°": np.tan(np.deg2rad(63.75)),
    }
    # Bow shock intercept
    inter_med = fitting_coeffs[6]**2 - (fitting_coeffs[0]*fitting_coeffs[9])
    x_s = (-fitting_coeffs[6] + np.sqrt(inter_med)) / fitting_coeffs[0]
    x_e = 30
    y_s = 0
    slope = line_slopes[angle_key]
    y_e = -x_e * slope
    # angle line parameters: (x_s, x_e, y_s, y_e)
    angle_line = (x_s, x_e, y_s, y_e)
    return(angle_line)


def draw_hist(ax, hist, extent, cmap, v_bounds, angle_line):
    """Draw heatmap + cone angle line. v_bounds is list taking lower, upper values in order."""
    ax.imshow(hist, interpolation='none', origin='lower',
              extent=extent, cmap=cmap, vmin=v_bounds[0], vmax=v_bounds[1])  

    x_s, x_e, y_s, y_e = angle_line
    ax.plot([x_s, x_e], [y_s, y_e], color='k', linewidth=1)

def draw_heatmap(ax, hist, extent, cmap, cmap_norm, angle_line):
    """Draw heatmap + flow line"""
    ax.imshow(hist, interpolation='none', origin='lower',
              extent=extent, cmap=cmap, norm=cmap_norm)
    
    x_s, x_e, y_s, y_e = angle_line
    ax.plot([x_s, x_e], [y_s, y_e], color='k', linewidth=1)

def set_limits(ax):
    ax.set_xlim(0, 20)
    ax.set_ylim(-20, 20)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.set_aspect('equal')

def mask_inside_magnetopause(ax, x_shue, r_shue):
    """
    Mask (fill white) the region inside the magnetopause and bounded by x=0.
    
    Parameters
    ----------
    ax : matplotlib axis
        Axis to draw on
    x_shue : array
        X coordinates of magnetopause
    r_shue : array
        R (Y) coordinates of magnetopause
    zorder : int zorder=10
        Draw order (should be higher than background contours)
    """

    # Ensure arrays are numpy arrays
    x_shue = np.asarray(x_shue)
    r_shue = np.asarray(r_shue)

    # Magnetopause runs from theta = -pi/2 to pi/2
    # So first point is lower flank, last point is upper flank
    x_lower, y_lower = x_shue[0], r_shue[0]
    x_upper, y_upper = x_shue[-1], r_shue[-1]

    # Build closed polygon:
    # magnetopause curve
    poly_x = list(x_shue)
    poly_y = list(r_shue)

    # connect upper flank to (0,0)
    poly_x.append(0)
    poly_y.append(0)

    # connect (0,0) back to lower flank
    poly_x.append(0)
    poly_y.append(0)

    # close polygon automatically by fill , zorder=zorder
    ax.fill(poly_x, poly_y, color='white')


######PRODUCING PLOTS!!!

def CA_MA_Binned_Plot(property_key,comparison_key, ma_ca_blocks, xedg, yedg):
    """Cone Angle and Mach Number 3x5 grid """
    #find constants for Merka & Shue
    X_shue, R_shue, Xgipm, Ygipm, Zgipm, f, fitting_coeffs = model_calcs()

    #main cbar titles:

    cbar_titles = {"Peak Transverse Frequency": "Frequency, mHz", "Peak Compressive Frequency": "Frequency, mHz", "Ellipticity": "Ellipticity", "Normalised Transverse Frequency": "Normalised Frequency, Hz/nT", "Normalised Compressive Frequency": "Normalised Frequency, Hz/nT"}

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

    cmap_dict = {"Peak Transverse Frequency": [0.007, 0.1], "Peak Compressive Frequency": [0.007, 0.1], "Ellipticity": [1,4], "Normalised Transverse Frequency": [0.001, 0.1], "Normalised Compressive Frequency": [0.001, 0.1],}
    
    if property_key=='Ellipticity':
        powercmp = cm_cram.navia
        power_norm = colors.Normalize(vmin=cmap_dict[property_key][0], vmax=cmap_dict[property_key][1])
    else:
        powercmp = cm_cram.lipari
        powercmp.set_under(powercmp(0))
        power_norm = colors.LogNorm(vmin=cmap_dict[property_key][0], vmax=cmap_dict[property_key][1])
    
    comparison_cmp = cm_cram.vik
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
    if property_key=='Ellipticity':
        cbar1.set_ticks(ticks=[1, 2, 3, 4], labels=['1', '2', '3', '4'])
    if property_key=="Peak Transverse Frequency":
        cbar1.set_ticks(ticks=[0.007, 0.01, 0.1], labels=['7', '10', '100'])
    if property_key=="Peak Compressive Frequency":
        cbar1.set_ticks(ticks=[0.007, 0.01, 0.1], labels=['7', '10', '100'])
    cbar2.ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext())

    path = "/Users/roseatkinson/Documents/New_Figs/CA_MA_" + property_key + comparison_key + ".png"
    plt.savefig(path)
    #Save Plots!!


##################
##################
##################
#Wider CA bin plot, with transverse & compressive in same batch

def NaNp_heatmap_plot(property_key, nanp_ca_blocks, xedg, yedg):
    """Plot 3 x 4 graph to show wave power"""
    #properties are either PEAK FREQUENCY or NORMALISED FREQUENCY
    #new angle titles!
    
    #find constants for Merka & Shue
    X_shue, R_shue, Xgipm, Ygipm, Zgipm, f, fitting_coeffs = model_calcs()

    angle_titles_wide = ["30-52.5°", "52.5-75°", "30-52.5°", "52.5-75°"]
    cbar_titles = {"Peak Frequency": "Frequency, mHz", "Normalised Frequency": "Normalised Frequency, Hz/nT"}
    extent = [xedg[0], xedg[-1], yedg[0], yedg[-1]]
    # -------------------------------
    # CREATE FIGURE + GRID
    # -------------------------------
    
    fig = plt.figure(figsize=(9, 8), dpi=300, constrained_layout=True)
    gs = fig.add_gridspec(
        nrows=3, ncols=5,      # 1 column for patch labels
        width_ratios=[0.35, 1, 1, 1, 1],  # label column thinner
        wspace=0.05, hspace=0.1
    )
    
    fig.suptitle("Effect of Alpha/Proton Ratio", fontsize=18)
    plt.rcParams['axes.labelsize'] = 14
    
    # Row labels (top row → bottom row)
    row_labels = [
        r'$0.035 \leq N_a/N_p$',
        r'$N_a/N_p < 0.025$',
        r'Ratio'
    ]
    # Line slopes for different angle classes
    line_slopes = {
        "0–30°": np.tan(np.deg2rad(15)),
        "30–45°": np.tan(np.deg2rad(37.5)),
        "45–60°": np.tan(np.deg2rad(52.5)),
        "60–75°": np.tan(np.deg2rad(67.5)),
        "75–90°": np.tan(np.deg2rad(82.5)),
        "30-52.5°": np.tan(np.deg2rad(41.25)),
        "52.5-75°": np.tan(np.deg2rad(63.75)),
    }

    # -------------------------------
    # MAKE AXES FOR THE 3×5 PANELS
    # -------------------------------
    axs = []
    for r in range(3):
        plot_row_axes = []
        gs_row = r      
        for c in range(4):
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

    cmap_dict = {"Peak Frequency": [0.007, 0.1], "Normalised Frequency": [0.001, 0.1]}

    powercmp = cm_cram.lipari
    #powercmp.set_under(powercmp(0))
    power_norm = colors.LogNorm(vmin=cmap_dict[property_key][0], vmax=cmap_dict[property_key][1])

    comparison_cmp = cm_cram.vik
    comparison_norm = colors.LogNorm(vmin=0.1, vmax=10)
    
    # -------------------------------
    # PLOT ALL PANELS
    # -------------------------------
    
    for col in range(4):                         # angle class
        title = angle_titles_wide[col]
        slope = line_slopes[title]
    
        for row in range(3):                     # mach no. class
            ax = axs[row][col]
    
            # Draw contour, magnetopause
            draw_background(ax, Xgipm[:, :, 0], Ygipm[:, :, 0], f[:, :, 0],
                            X_shue, R_shue)
    
            # Histogram for this cell
            hist = nanp_ca_blocks[row][col]
    
            # angle line parameters: (x_s, x_e, y_s, y_e)
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
    
    # --- Top two rows colourbar (wave power)
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
    cbar2.set_label('Ratio')
    
    cbar1.ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext())
    cbar2.ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext())

    path = "/Users/roseatkinson/Documents/New_Figs/CA_NaNp" + property_key +".png"
    plt.savefig(path)

##################
##################
##################
#Wider CA bin plot, with transverse & compressive in same batch

def NaNp_heatmap_plot(property_key, nanp_ca_blocks, xedg, yedg, NaNp_lims):
    """Plot 3 x 4 graph to show wave power"""
    #properties are either PEAK FREQUENCY or NORMALISED FREQUENCY
    #new angle titles!
    
    #find constants for Merka & Shue
    X_shue, R_shue, Xgipm, Ygipm, Zgipm, f, fitting_coeffs = model_calcs()

    angle_titles_wide = ["30-52.5°", "52.5-75°", "30-52.5°", "52.5-75°"]
    cbar_titles = {"Peak Frequency": "Frequency, mHz", "Normalised Frequency": "Normalised Frequency, Hz/nT"}
    extent = [xedg[0], xedg[-1], yedg[0], yedg[-1]]
    # -------------------------------
    # CREATE FIGURE + GRID
    # -------------------------------
    
    fig = plt.figure(figsize=(9, 8), dpi=300, constrained_layout=True)
    gs = fig.add_gridspec(
        nrows=3, ncols=5,      # 1 column for patch labels
        width_ratios=[0.35, 1, 1, 1, 1],  # label column thinner
        wspace=0.05, hspace=0.1
    )
    
    fig.suptitle("Effect of Alpha/Proton Ratio", fontsize=18)
    plt.rcParams['axes.labelsize'] = 14
    
    # Row labels (top row → bottom row)
    row_labels = [
        NaNp_lims[1] + r'$ \leq N_a/N_p$',
        r'$N_a/N_p < $' + NaNp_lims[0],
        r'Ratio'
    ]
    # Line slopes for different angle classes
    line_slopes = {
        "0–30°": np.tan(np.deg2rad(15)),
        "30–45°": np.tan(np.deg2rad(37.5)),
        "45–60°": np.tan(np.deg2rad(52.5)),
        "60–75°": np.tan(np.deg2rad(67.5)),
        "75–90°": np.tan(np.deg2rad(82.5)),
        "30-52.5°": np.tan(np.deg2rad(41.25)),
        "52.5-75°": np.tan(np.deg2rad(63.75)),
    }

    # -------------------------------
    # MAKE AXES FOR THE 3×5 PANELS
    # -------------------------------
    axs = []
    for r in range(3):
        plot_row_axes = []
        gs_row = r      
        for c in range(4):
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

    cmap_dict = {"Peak Frequency": [0.007, 0.1], "Normalised Frequency": [0.001, 0.1]}

    powercmp = cm_cram.lipari
    #powercmp.set_under(powercmp(0))
    power_norm = colors.LogNorm(vmin=cmap_dict[property_key][0], vmax=cmap_dict[property_key][1])

    comparison_cmp = cm_cram.vik
    comparison_norm = colors.LogNorm(vmin=0.1, vmax=10)
    
    # -------------------------------
    # PLOT ALL PANELS
    # -------------------------------
    
    for col in range(4):                         # angle class
        title = angle_titles_wide[col]
        slope = line_slopes[title]
    
        for row in range(3):                     # mach no. class
            ax = axs[row][col]
    
            # Draw contour, magnetopause
            draw_background(ax, Xgipm[:, :, 0], Ygipm[:, :, 0], f[:, :, 0],
                            X_shue, R_shue)
    
            # Histogram for this cell
            hist = nanp_ca_blocks[row][col]
    
            # angle line parameters: (x_s, x_e, y_s, y_e)
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
    
    # --- Top two rows colourbar (wave power)
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
    cbar2.set_label('Ratio')
    
    cbar1.ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext())
    cbar2.ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext())

    path = "/Users/roseatkinson/Documents/New_Figs/CA_NaNp" + property_key + "lims_" + NaNp_lims[0] + NaNp_lims[1] +".png"
    plt.savefig(path)

def NaNp_heatmap_plot_B_filt(property_key, nanp_ca_blocks, xedg, yedg, NaNp_lims):
    """Plot 3 x 4 graph to show wave power"""
    #properties are either PEAK FREQUENCY or NORMALISED FREQUENCY
    #new angle titles!
    
    #find constants for Merka & Shue
    X_shue, R_shue, Xgipm, Ygipm, Zgipm, f, fitting_coeffs = model_calcs()

    angle_titles_wide = ["30-52.5°", "52.5-75°", "30-52.5°", "52.5-75°"]
    cbar_titles = {"Peak Frequency": "Frequency, mHz", "Normalised Frequency": "Normalised Frequency, Hz/nT"}
    extent = [xedg[0], xedg[-1], yedg[0], yedg[-1]]
    # -------------------------------
    # CREATE FIGURE + GRID
    # -------------------------------
    
    fig = plt.figure(figsize=(9, 8), dpi=300, constrained_layout=True)
    gs = fig.add_gridspec(
        nrows=3, ncols=5,      # 1 column for patch labels
        width_ratios=[0.35, 1, 1, 1, 1],  # label column thinner
        wspace=0.05, hspace=0.1
    )
    
    fig.suptitle("Effect of Alpha/Proton Ratio, Constant B", fontsize=18, y=1.08)
    plt.rcParams['axes.labelsize'] = 14
    
    # Row labels (top row → bottom row)
    row_labels = [
        NaNp_lims[1] + r'$ \leq N_a/N_p$',
        r'$N_a/N_p < $' + NaNp_lims[0],
        r'Ratio'
    ]
    # Line slopes for different angle classes
    line_slopes = {
        "0–30°": np.tan(np.deg2rad(15)),
        "30–45°": np.tan(np.deg2rad(37.5)),
        "45–60°": np.tan(np.deg2rad(52.5)),
        "60–75°": np.tan(np.deg2rad(67.5)),
        "75–90°": np.tan(np.deg2rad(82.5)),
        "30-52.5°": np.tan(np.deg2rad(41.25)),
        "52.5-75°": np.tan(np.deg2rad(63.75)),
    }

    # -------------------------------
    # MAKE AXES FOR THE 3×5 PANELS
    # -------------------------------
    axs = []
    for r in range(3):
        plot_row_axes = []
        gs_row = r      
        for c in range(4):
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

    cmap_dict = {"Peak Frequency": [0.007, 0.1], "Normalised Frequency": [0.001, 0.1]}

    powercmp = cm_cram.lipari
    #powercmp.set_under(powercmp(0))
    power_norm = colors.LogNorm(vmin=cmap_dict[property_key][0], vmax=cmap_dict[property_key][1])

    comparison_cmp = cm_cram.vik
    comparison_norm = colors.LogNorm(vmin=0.1, vmax=10)
    
    # -------------------------------
    # PLOT ALL PANELS
    # -------------------------------
    
    for col in range(4):                         # angle class
        title = angle_titles_wide[col]
        slope = line_slopes[title]
    
        for row in range(3):                     # mach no. class
            ax = axs[row][col]
    
            # Draw contour, magnetopause
            draw_background(ax, Xgipm[:, :, 0], Ygipm[:, :, 0], f[:, :, 0],
                            X_shue, R_shue)
    
            # Histogram for this cell
            hist = nanp_ca_blocks[row][col]
    
            # angle line parameters: (x_s, x_e, y_s, y_e)
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
    # COLUMN GROUP LABELS
    # -------------------------------
    
    # Ensure layout is computed so positions are correct
    fig.canvas.draw()
    
    # Get positions of top row axes
    pos0 = axs[0][0].get_position()
    pos1 = axs[0][1].get_position()
    pos2 = axs[0][2].get_position()
    pos3 = axs[0][3].get_position()
    
    # Midpoints for the two column groups
    x_transverse = ((pos0.x0 + pos1.x1) / 2)-0.02
    x_compressive = ((pos2.x0 + pos3.x1) / 2)-0.08
    
    # Vertical position slightly above column titles
    y_top = pos0.y1 + 0.035
    
    # Add figure-level text
    fig.text(x_transverse, y_top, "Transverse",
             ha='center', va='bottom', fontsize=14, fontweight='bold')
    
    fig.text(x_compressive, y_top, "Compressive",
             ha='center', va='bottom', fontsize=14, fontweight='bold')
    
    # -------------------------------
    # COLORBARS (TWO SEPARATE ON RIGHT)
    # -------------------------------
    
    from matplotlib.cm import ScalarMappable
    
    # --- Scalar mappables (independent of any single subplot image)
    sm_power = ScalarMappable(norm=power_norm, cmap=powercmp)
    sm_power.set_array([])
    
    sm_ratio = ScalarMappable(norm=comparison_norm, cmap=comparison_cmp)
    sm_ratio.set_array([])
    
    # --- Top two rows colourbar (wave power)
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
    cbar2.set_label('Ratio')
    
    cbar1.ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext())
    cbar2.ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext())

    path = "/Users/roseatkinson/Documents/New_Figs/CA_NaNp" + property_key + "lims_" + NaNp_lims[0] + NaNp_lims[1] + "B_const_3_5.png"
    plt.savefig(path)


############################


def CA_Error_Plot(property_key, ca_blocks, xedg, yedg):
    """Cone Angle 1x5 grid """
    #find constants for Merka & Shue
    X_shue, R_shue, Xgipm, Ygipm, Zgipm, f, fitting_coeffs = model_calcs()

    fig = plt.figure(figsize=(8.5, 3), dpi=300, constrained_layout=True)
    gs = fig.add_gridspec(
        nrows=1, ncols=6,      # 1 column for patch labels
        width_ratios=[0.35, 1, 1, 1, 1, 1],  # label column thinner
        wspace=0.05, hspace=0.1
    )

    fig.suptitle(property_key, fontsize=18)
    plt.rcParams['axes.labelsize'] = 14
    
    angle_titles = ["0–30°", "30–45°", "45–60°", "60–75°", "75–90°"]
    extent = [xedg[0], xedg[-1], yedg[0], yedg[-1]]
    # -------------------------------
    # MAKE AXES FOR THE 1×5 PANELS
    # -------------------------------
    
    axs = []     
    for c in range(5):
        ax = fig.add_subplot(gs[c + 1])
        axs.append(ax)

    # -------------------------------
    # COLORMAP
    # -------------------------------

    powercmp = cm_cram.vik
    power_norm = colors.TwoSlopeNorm(vmin=-5, vcenter=0., vmax=5)

    # -------------------------------
    # PLOT ALL PANELS
    # -------------------------------

    for col in range(5):    # angle class
        
        title = angle_titles[col]
        ax = axs[col]

        # Draw contour, magnetopause
        draw_background(ax, Xgipm[:, :, 0], Ygipm[:, :, 0], f[:, :, 0],
                        X_shue, R_shue)

        # Histogram for this cell
        hist = ca_blocks[col]
        
        angle_line = cone_angle_line(fitting_coeffs, title)

        draw_heatmap(ax, hist, extent, powercmp, power_norm, angle_line)
            
        mask_inside_magnetopause(ax, X_shue, R_shue)
            
        # redraw magnetopause boundary so it stays crisp
        ax.plot(X_shue, R_shue, 'k', linewidth=1)
        set_limits(ax)

        # Labels
        ax.set_title(rf'$\alpha$ = {title}', fontsize=12)
        ax.set_xlabel("$X_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")
        if col == 0:
            ax.set_ylabel("$Y_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")
            
            

    # -------------------------------
    # COLORBAR
    # -------------------------------

    from matplotlib.cm import ScalarMappable

    # --- Scalar mappables (independent of any single subplot image)
    sm_power = ScalarMappable(norm=power_norm, cmap=powercmp)
    sm_power.set_array([])

    cbar1 = fig.colorbar(
        sm_power,
        ax=axs[4],
        location='right',
        pad=0.02,
        extend='both'
    )
    cbar1.set_label('Error')

    path = "/Users/roseatkinson/Documents/New_Figs/CA_" + property_key + ".png"
    plt.savefig(path)


###############################################################################

def CA_Abs_Error_Plot(property_key, ca_blocks, hist_blocks, xedg, yedg):
    """Cone Angle 1x5 grid """
    
    #find constants for Merka & Shue
    X_shue, R_shue, Xgipm, Ygipm, Zgipm, f, fitting_coeffs = model_calcs()

    fig = plt.figure(figsize=(8.5, 3), dpi=300, constrained_layout=True)
    gs = fig.add_gridspec(
        nrows=1, ncols=6,      # 1 column for patch labels
        width_ratios=[0.35, 1, 1, 1, 1, 1],  # label column thinner
        wspace=0.05, hspace=0.1
    )

    fig.suptitle(property_key, fontsize=18)
    plt.rcParams['axes.labelsize'] = 14
    
    angle_titles = ["0–30°", "30–45°", "45–60°", "60–75°", "75–90°"]
    extent = [xedg[0], xedg[-1], yedg[0], yedg[-1]]
    # -------------------------------
    # MAKE AXES FOR THE 1×5 PANELS
    # -------------------------------
    
    axs = []     
    for c in range(5):
        ax = fig.add_subplot(gs[c + 1])
        axs.append(ax)

    # -------------------------------
    # COLORMAP
    # -------------------------------
    
    #cmap = (matplotlib.colors.ListedColormap(['darkblue', 'lightblue', 'white', 'palevioletred', 'brown']).with_extremes(under='darkblue', over='brown'))
    cmap = matplotlib.colormaps['RdBu_r'].resampled(5)
    bounds = [-3, -2, -1, 1, 2, 3]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    # cmap_obs = (matplotlib.colors.ListedColormap(['lightgrey','lightgrey', 'lightgrey']))
    # bounds_obs = [0,10,20,50]
    # norm_obs = matplotlib.colors.BoundaryNorm(bounds_obs, cmap_obs.N)

    cmap_obs = (matplotlib.colors.ListedColormap(['black', 'grey', 'white', 'palevioletred', 'brown']).with_extremes(under='darkblue', over='brown'))
    bounds_obs = [0, 1, 2, 3, 4, 5]
    norm_obs = matplotlib.colors.BoundaryNorm(bounds_obs, cmap_obs.N)

    #plt.rc('image', cmap=cmap_obs)
    # -------------------------------
    # PLOT ALL PANELS
    # -------------------------------

    for col in range(5):    # angle class
        
        title = angle_titles[col]
        ax = axs[col]
        ax.fill([0 , 20, 20, 0], [20, 20, -20, -20], "lightgrey", zorder=0)
        # Draw contour, magnetopause
        draw_background(ax, Xgipm[:, :, 0], Ygipm[:, :, 0], f[:, :, 0],
                        X_shue, R_shue)

        # Histogram for this cell
        hist = ca_blocks[col]
        obs = hist_blocks[col]
        
        angle_line = cone_angle_line(fitting_coeffs, title)
        ax.imshow(obs, interpolation='none', origin='lower',
               extent=extent, cmap=cmap_obs, norm=norm_obs)
        draw_heatmap(ax, hist, extent, cmap, norm, angle_line)
            
        mask_inside_magnetopause(ax, X_shue, R_shue)
            
        # redraw magnetopause boundary so it stays crisp
        ax.plot(X_shue, R_shue, 'k', linewidth=1)
        set_limits(ax)

        # Labels
        ax.set_title(rf'$\alpha$ = {title}', fontsize=12)
        ax.set_xlabel("$X_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")
        if col == 0:
            ax.set_ylabel("$Y_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")
            

    # -------------------------------
    # COLORBAR
    # -------------------------------

    from matplotlib.cm import ScalarMappable

    # --- Scalar mappables (independent of any single subplot image)
    sm_power = ScalarMappable(cmap=cmap, norm=norm)
    sm_power.set_array([])

    cbar1 = fig.colorbar(
        sm_power,
        ax=axs[4],
        location='right',
        pad=0.02,
        spacing='proportional',
        extend='both'
    )
    
    cbar1.set_label('Error/Resolution')

    path = "/Users/roseatkinson/Documents/New_Figs/CA_" + property_key + ".png"
    plt.savefig(path)



#location map

def locationmapping(location_list, location_refs):
    """Unpack all points in location list as x,y lists and plot"""
    X_shue, R_shue, Xgipm, Ygipm, Zgipm, f, fitting_coeffs = model_calcs()

    ###################
    fig = plt.figure(figsize=(6, 4.5))
    subfigs = fig.subfigures(1, 1)
    axsLeft = subfigs.subplots()

    loc_x = []
    loc_y = []
    loc_label_list_x = []
    
    for loc_pair in location_list:
        loc_label_x = loc_pair[0]-1
        loc_label_list_x.append(loc_label_x)
        loc_x.append(loc_pair[0])
        loc_y.append(loc_pair[1])

    draw_background(axsLeft, Xgipm[:, :, 0], Ygipm[:, :, 0], f[:, :, 0],
                        X_shue, R_shue)
    
    angle_line = cone_angle_line(fitting_coeffs, "0–30°")
    
    x_s, x_e, y_s, y_e = angle_line
    
    axsLeft.plot([x_s, x_e], [y_s, y_e], color='k', linewidth=1)
    
    axsLeft.scatter(loc_x, loc_y, marker="o", c="blue", s=60)

    for loc_lab_x, loc_lab_y, loc_ref in zip(loc_label_list_x, loc_y, location_refs):
        plt.text(loc_lab_x, loc_lab_y, loc_ref, color="maroon", weight= 'bold', fontsize=16)

    set_limits(axsLeft)

    axsLeft.set_title(r'Case Study Locations')
    axsLeft.set_xlabel("$X_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")
    axsLeft.set_ylabel("$Y_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")

    path = "/Users/roseatkinson/Documents/New_Figs/AlphaCaseStudies.png"
    plt.savefig(path)