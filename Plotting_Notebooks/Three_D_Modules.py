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

    #m_1 = 10.22
    #m_2 = 1.29
    #m_3 = 0.184
    #m_4 = 8.14
    #m_5 = 6.6
    #m_6 = 0.58
    #m_7 = -0.007
    #m_8 = 0.024

    #use B_z = 0
    #B_z = 0
    #D_p = 1.76

    #alpha = (m_6 + m_7*B_z)*(1 +m_8*(np.log(D_p)))

    #tanh_angle = m_3*(B_z+m_4)
    #r_0 = (m_1 + m_2*np.tanh(tanh_angle))*(D_p**(-1/m_5))

    #pi = np.pi

    #theta = np.arange(-pi/2, pi/2, 0.01)

    #r_mod = (2/(1+np.cos(theta)))**alpha
    #r = r_0*r_mod

    #X_shue = r*(np.cos(theta))
    #R_shue = r*(np.sin(theta))

    x = np.linspace(0, 20, 100) #x coordinates (Re)
    y = np.linspace(-30, 30, 100) #y coordinates (Re)
    z_0 = 0 #z coordinates in Re
    z_1 = 2 #z coordinates in Re
    z_2 = 4 #z coordinates in Re
    z_3 = 6 #z coordinates in Re
    
    [Xgipm_0,Ygipm_0,Zgipm_0] = np.meshgrid(x,y,z_0,indexing="ij")
    [Xgipm_1,Ygipm_1,Zgipm_1] = np.meshgrid(x,y,z_1,indexing="ij")
    [Xgipm_2,Ygipm_2,Zgipm_2] = np.meshgrid(x,y,z_2,indexing="ij")
    [Xgipm_3,Ygipm_3,Zgipm_3] = np.meshgrid(x,y,z_3,indexing="ij")
    
    XMA_all = 10
    fitting_coeffs = merka05_surface_eq_array_GIPM(XMA_all)

    Xn_0 = Xgipm_0
    Yn_0 = Ygipm_0
    Zn_0 = Zgipm_0
    f_0 = fitting_coeffs[0]*Xn_0**2 + fitting_coeffs[1]*Yn_0**2 + fitting_coeffs[2]*Zn_0**2+ 2*fitting_coeffs[3]*Xn_0*Yn_0 + 2*fitting_coeffs[4]*Yn_0*Zn_0 + 2*fitting_coeffs[5]*Xn_0*Zn_0 + 2*fitting_coeffs[6]*Xn_0+2*fitting_coeffs[7]*Yn_0 + 2*fitting_coeffs[8]*Zn_0 + fitting_coeffs[9]

    set_0 = [Xn_0, Yn_0, Zn_0, f_0]

    Xn_1 = Xgipm_1
    Yn_1 = Ygipm_1
    Zn_1 = Zgipm_1
    f_1 = fitting_coeffs[0]*Xn_1**2 + fitting_coeffs[1]*Yn_1**2 + fitting_coeffs[2]*Zn_1**2+ 2*fitting_coeffs[3]*Xn_1*Yn_1 + 2*fitting_coeffs[4]*Yn_1*Zn_1 + 2*fitting_coeffs[5]*Xn_1*Zn_1 + 2*fitting_coeffs[6]*Xn_1+2*fitting_coeffs[7]*Yn_1 + 2*fitting_coeffs[8]*Zn_1 + fitting_coeffs[9]

    set_1 = [Xn_1, Yn_1, Zn_1, f_1]
    
    Xn_2 = Xgipm_2
    Yn_2 = Ygipm_2
    Zn_2 = Zgipm_2
    f_2 = fitting_coeffs[0]*Xn_2**2 + fitting_coeffs[1]*Yn_2**2 + fitting_coeffs[2]*Zn_2**2+ 2*fitting_coeffs[3]*Xn_2*Yn_2 + 2*fitting_coeffs[4]*Yn_2*Zn_2 + 2*fitting_coeffs[5]*Xn_2*Zn_2 + 2*fitting_coeffs[6]*Xn_2+2*fitting_coeffs[7]*Yn_2 + 2*fitting_coeffs[8]*Zn_2 + fitting_coeffs[9]

    set_2 = [Xn_2, Yn_2, Zn_2, f_2]
    
    Xn_3 = Xgipm_3
    Yn_3 = Ygipm_3
    Zn_3 = Zgipm_3
    f_3 = fitting_coeffs[0]*Xn_3**2 + fitting_coeffs[1]*Yn_3**2 + fitting_coeffs[2]*Zn_3**2+ 2*fitting_coeffs[3]*Xn_3*Yn_3 + 2*fitting_coeffs[4]*Yn_3*Zn_3 + 2*fitting_coeffs[5]*Xn_3*Zn_3 + 2*fitting_coeffs[6]*Xn_3+2*fitting_coeffs[7]*Yn_3 + 2*fitting_coeffs[8]*Zn_3 + fitting_coeffs[9]

    set_3 = [Xn_3, Yn_3, Zn_3, f_3]
    
    return(set_0, set_1, set_2, set_3, fitting_coeffs)


def draw_background(ax, coord_set, level):
    """Draw bow shock, magnetopause, and y=0 line"""
    ax.contour(coord_set[0], coord_set[1], coord_set[3], levels=[level], colors="black", linewidths=1)
    #ax.plot(x_shue, r_shue, 'k', linewidth=1)
    ax.hlines(0, 0, 25, color='k', linewidth=1)

def cone_angle_line(fitting_coeffs, angle_key):
    # Line slopes for different angle classes
    line_slopes = {
        "0–30°": np.tan(np.deg2rad(15)),
        "20–40°": np.tan(np.deg2rad(30)),
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

def ThreeD_Binned_Plot(property_key, threeD_blocks, xedg, yedg):
    """4-step Z 0 to 6 grid"""
    #find constants for Merka & Shue
    level_0_set, level_1_set, level_2_set, level_3_set, fitting_coeffs = model_calcs()

    #main cbar titles:

    cbar_titles = {"Transverse Power": "Normalised Wave Power", "Peak Compressive Frequency": "Normalised Wave Power"}

    fig = plt.figure(figsize=(9, 8), dpi=300, constrained_layout=True)
    gs = fig.add_gridspec(
        nrows=4, ncols=2,      # 1 column for patch labels
        width_ratios=[0.35, 1],  # label column thinner
        wspace=0.05, hspace=0.1
    )

    fig.suptitle(property_key, fontsize=18)
    plt.rcParams['axes.labelsize'] = 14

    # Row labels (top row → bottom row)
    row_labels = [
        r'$Z_{GIPM} = 6 R_E$',
        r'$Z_{GIPM} = 4 R_E$',
        r'$Z_{GIPM} = 2 R_E$',
        r'$Z_{GIPM} = 0 R_E$',
    ]
    extent = [xedg[0], xedg[-1], yedg[0], yedg[-1]]
    # -------------------------------
    # MAKE AXES FOR THE PANELS
    # -------------------------------
    
    axs = []
    for r in range(4):
        gs_row = r      
        ax = fig.add_subplot(gs[gs_row, 1])
        axs.append(ax)

    # -------------------------------
    # Patch Labels (Rounded Boxes)
    # -------------------------------

    for r in range(4):
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

    powercmp = 'magma'
    power_norm = colors.LogNorm(vmin=0.001, vmax=0.1)

    # -------------------------------
    # PLOT ALL PANELS
    # -------------------------------

    for row in range(4):                     # mach no. class
        ax = axs[row]

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
            
        #mask_inside_magnetopause(ax, X_shue, R_shue)
        
        # redraw magnetopause boundary so it stays crisp
        ax.plot(X_shue, R_shue, 'k', linewidth=1)
        set_limits(ax)

        # Labels
        ax.set_ylabel("$Y_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")
        if row == 0:
            ax.set_title(rf'$\alpha$ = {title}', fontsize=12)
        if row == 4:
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

    # --- Colourbar
    top_axes = axs[1] + axs[2]   # flatten row 0 and 1
    cbar1 = fig.colorbar(
        sm_power,
        ax=top_axes,
        location='right',
        pad=0.02,
        extend='both'
    )
    cbar1.set_label(cbar_titles[property_key])

    cbar1.ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext())

    path = "/Users/roseatkinson/Documents/New_Figs/3D_Plots_" + property_key + ".png"
    plt.savefig(path)
    #Save Plots!!

