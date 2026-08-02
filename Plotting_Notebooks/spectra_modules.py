import pandas as pd
import numpy as np
import datetime as dt

import matplotlib 
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

def spectra_plot(f_df, dt_str, freq_dict, alpha_val, sc):
    
    frequency = f_df['Frequency']
    perp_power = f_df['Total Transverse Power']
    para_power = f_df['Compressive Power']
    
    int_lower_lim = 7*(10**(-3))
    int_upper_lim = 100*(10**(-3))

    fig, ax = plt.subplots(figsize = (5, 8))
    
    ax.set_title(r'$N_{\mathrm{\alpha}} / N_{\mathrm{p}} =$' + str(alpha_val))
    ax.plot(frequency, perp_power, color="black", label='Compressive Power')
    ax.plot(frequency, para_power, color="red", label='Transverse Power')
    #ax.plot(x_array, y_array, color="blue", label='-5/3 Slope')
    #plt.xlabel('Period (s)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'Power, $\mathrm{nT^2/Hz}$')
    ax.set_xlabel('Frequency, Hz')
    ax.set_xlim(0.001, 5)
    ax.set_ylim(0.000001, 1000)
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
    #ax.yaxis.set_minor_locator(ticker.LogLocator(base=10, numticks=10))
    ax.vlines(x=[int_lower_lim, int_upper_lim], ymin = 0.0000001, ymax = 10000, linestyles='dotted', color='mediumblue', label='ULF Band')
    ax.vlines(x=freq_dict['peak_freq'], ymin = 0.0000001, ymax = 100_000_000, linestyles='solid', color='mediumblue', label='Peak Frequency')
    ax.vlines(x=freq_dict['tak_freq'], ymin = 0.0000001, ymax = 100_000_000, linestyles='dashed', color='mediumblue', label='Takahashi Frequency')
    #ax.vlines(x=int_upper_lim, ymin = 0.0000001, ymax = 10000, linestyles='dotted', color='mediumblue', label='ULF Upper Bound')
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=14) 
    plt.rc('axes', labelsize=16) 
    
    ax.set_axisbelow(True)
    ax.yaxis.grid(color='lightgray')
    ax.xaxis.grid(color='lightgray')
    
    path = "/Users/roseatkinson/Documents/New_Figs/Nanp"+str(alpha_val)+sc+dt_str+".png"
    plt.savefig(path, bbox_inches='tight')

################################################################################

def TS_plot(csv_path , str_centre):

    time_cent = pd.to_datetime(str_centre)

    df = pd.read_csv(csv_path)
    df['datetime'] = pd.to_datetime(df['datetime'])
    df = df.set_index('datetime')

    locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
    formatter = mdates.ConciseDateFormatter(locator)

    #mask to just relevant times
    time_start_10mins = time_cent - dt.timedelta(minutes=5)
    time_end_10mins = time_cent + dt.timedelta(minutes=5)
    
    time_start_4mins = time_cent - dt.timedelta(minutes=2)
    time_end_4mins = time_cent + dt.timedelta(minutes=2)

    df = df.loc[((df.index >= time_start_10mins) & (df.index < time_end_10mins))]
    
    fig, ax = plt.subplots(figsize = (5, 8))
    
    ax.set_title('Original Time Series')
    ax.plot(df.index, df['Bx_gse'], color = 'g', label='$B_{x}$')
    ax.plot(df.index, df['By_gse'], color = 'b', label='$B_{y}$')
    ax.plot(df.index, df['Bz_gse'], color = 'r', label='$B_{z}$')
    ax.plot(df.index, df['B_mag'], color = 'k', label='$B_{mag}$')
    ax.vlines(x=[time_start_4mins, time_end_4mins], ymin = -10, ymax = 10, linestyles='dotted', color='k', label='ULF Interval')
    ax.set_xlabel('Time')
    ax.set_ylabel('Magnetic Field, nT')
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    ax.set_axisbelow(True)
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.grid(color='lightgray')
    ax.xaxis.grid(color='lightgray')


#########################################################################

def TS_and_FS(ts_dfs, fs_dfs, labels, **kwargs):

    #ts_dfs, fs_dfs, labels should be a list, ts dfs should be shortened
    power_dict = {'abs': r'Power Spectral Density, $\mathrm{nT^2/Hz}$', 'norm': r'Power Spectral Density, $\mathrm{nT^2/Hz}$'}

    if 'norm' in kwargs:
        if kwargs['norm']==True:
            y_label = power_dict['norm']
        else:
            y_label = power_dict['abs']
    else:
        y_label = power_dict['abs']

    fig_width = 4*len(ts_dfs)
    
    # Create a figure
    fig = plt.figure(figsize=(fig_width, 9))
    
    locator = mdates.AutoDateLocator(minticks=3, maxticks=5)
    formatter = mdates.ConciseDateFormatter(locator)
    
    # Define a gridspec layout
    gs = GridSpec(4, len(ts_dfs), figure=fig)

    #Define ULF wave band

    int_lower_lim = 7*(10**(-3))
    int_upper_lim = 100*(10**(-3))
    
    # Add subplots with custom spans
    
    #Subplots are positioned and sized using slicing 
    #(e.g., gs[0, 0:2] spans the first two columns of the first row).

    ts_axes = []
    fs_axes = []

    for c in range(len(ts_dfs)):
        col = c
        ts_ax = fig.add_subplot(gs[0, c])
        fs_ax = fig.add_subplot(gs[1:, c])
        ts_axes.append(ts_ax)
        fs_axes.append(fs_ax)

    for df,ax, title in zip(ts_dfs,ts_axes, labels):
        
        # Add data to each subplot
        ax.plot(df.index, df['Bx_gse'], color = 'red')
        ax.plot(df.index, df['By_gse'], color = 'green')
        ax.plot(df.index, df['Bz_gse'], color = 'blue')
        ax.plot(df.index, df['B_mag'], color = 'k')
        ax.set_ylim(-10, 10)
        ax.set_xlim(df.index[0], df.index[-1])
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.xaxis.set_tick_params(labelsize=12)
        ax.yaxis.set_tick_params(labelsize=12)
        ax.set_ylabel(r'B, nT')
        ax.set_title(title)
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='lightgray')
        ax.xaxis.grid(color='lightgray')

    if 'tak_freq' in kwargs:
                
        for df,ax,freq in zip(fs_dfs,fs_axes, kwargs['tak_freq']):

            ax.plot(df['Frequency'], df['Total Transverse Power'], color="black", label='Transverse')
            ax.plot(df['Frequency'], df['Compressive Power'], color="red", label='Compressive')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_ylabel(y_label)
            ax.set_xlabel('Frequency, Hz')
            ax.set_xlim(0.001, 5)
            ax.set_ylim(0.000001, 1000)
            ax.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
            ax.vlines(x=int_lower_lim, ymin = 0.0000001, ymax = 10000, linestyles='dotted', color='k')
            ax.vlines(x=int_upper_lim, ymin = 0.0000001, ymax = 10000, linestyles='dotted', color='k')
            ax.axvspan(int_lower_lim, int_upper_lim, color='orange', alpha=0.5)
            ax.vlines(x=freq, ymin = 0.0000001, ymax = 100_000_000, linestyles='solid', color='k', label='Proton (Takahashi)')
            ax.vlines(x=freq/2, ymin = 0.0000001, ymax = 100_000_000, linestyles='dashed', color='k', label='Alpha (Takahashi)')
            ax.xaxis.set_tick_params(labelsize=12)
            ax.yaxis.set_tick_params(labelsize=12)
            ax.set_axisbelow(True)
            ax.yaxis.grid(color='lightgray')
            ax.xaxis.grid(color='lightgray')
            ax.legend(loc='lower left')

    else:
        for df,ax in zip(fs_dfs,fs_axes):
            ax.plot(df['Frequency'], df['Total Transverse Power'], color="black", label='Transverse')
            ax.plot(df['Frequency'], df['Compressive Power'], color="red", label='Compressive')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_ylabel(y_label)
            ax.set_xlabel('Frequency, Hz')
            ax.set_xlim(0.001, 5)
            ax.set_ylim(0.000001, 1000)
            ax.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
            ax.vlines(x=int_lower_lim, ymin = 0.0000001, ymax = 10000, linestyles='dotted', color='k')
            ax.vlines(x=int_upper_lim, ymin = 0.0000001, ymax = 10000, linestyles='dotted', color='k')
            ax.axvspan(int_lower_lim, int_upper_lim, color='orange', alpha=0.5)
            ax.xaxis.set_tick_params(labelsize=12)
            ax.yaxis.set_tick_params(labelsize=12)
            ax.set_axisbelow(True)
            ax.yaxis.grid(color='lightgray')
            ax.xaxis.grid(color='lightgray')
            ax.legend(loc='upper right')

    plt.rcParams['axes.titlesize'] = 18
    plt.rcParams['axes.labelsize'] = 18
    
    # Adjust layout
    plt.tight_layout()

    if 'filename' in kwargs:
        fname = kwargs['filename']
        path = "/Users/roseatkinson/Documents/New_Figs/Nanp"+fname+".png"
    else:
        path = "/Users/roseatkinson/Documents/New_Figs/Nanp_Temp.png"
        
    plt.savefig(path, bbox_inches='tight')
    
#############################################################

def spectra_overlay_plot(f_df_list, alpha_bin, cone_angle):
    
    int_lower_lim = 7*(10**(-3))
    int_upper_lim = 100*(10**(-3))

    fig, ax = plt.subplots(figsize = (5, 8))
    
    ax.set_title(r"$N_{\mathrm{\alpha}} / N_{\mathrm{p}}$" + ' ' + alpha_bin + ', cone angle ' + cone_angle)
    
    for f_df in f_df_list:
        ax.plot(f_df['Frequency'], f_df['Total Transverse Power'], color="dimgrey")
        ax.plot(f_df['Frequency'], f_df['Compressive Power'], color="palevioletred")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'Power, $\mathrm{nT^2/Hz}$')
    ax.set_xlabel('Frequency, Hz')
    ax.set_xlim(0.001, 5)
    ax.set_ylim(0.000001, 1000)
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
    #ax.yaxis.set_minor_locator(ticker.LogLocator(base=10, numticks=10))
    ax.vlines(x=[int_lower_lim, int_upper_lim], ymin = 0.0000001, ymax = 10000, linestyles='dotted', color='black', label='ULF Band')
    #ax.vlines(x=int_upper_lim, ymin = 0.0000001, ymax = 10000, linestyles='dotted', color='mediumblue', label='ULF Upper Bound')
    colors = ['dimgrey', 'palevioletred']
    lines = [Line2D([0], [0], color=c, linewidth=2) for c in colors]
    labels = ['Transverse', 'Compressive']
    plt.legend(lines, labels)
    #ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=14) 
    plt.rc('axes', labelsize=16) 
    
    ax.set_axisbelow(True)
    ax.yaxis.grid(color='lightgray')
    ax.xaxis.grid(color='lightgray')

    path = "/Users/roseatkinson/Documents/New_Figs/NanpComp"+str(alpha_bin)+cone_angle+".png"
    plt.savefig(path, bbox_inches='tight')
