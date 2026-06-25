import pandas as pd
import numpy as np
import datetime as dt

import matplotlib 
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as ticker

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
    plt.savefig(path)