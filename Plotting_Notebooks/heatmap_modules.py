import numpy as np 
import pandas as pd
from scipy import stats

##################################################################

def compute_hists2d(df, metric, *keys, **choices):
    """Compute 2D heatmaps of transverse and compressive power/frequency/compressibilty as specified, ignoring bins with obs_min no. of intervals (50 default)"""

    x_bin_edges = range(20)
    y_bin_edges = range(-20, 20)
    
    x_col='GIPM X (OMNI mean)'
    y_col='GIPM Y (OMNI mean)'

    w_dict = {'Transverse Power':'ULF Band Normalised Transverse Power', 'Compressive Power':'ULF Band Normalised Compressive Power', 'Compressibility':'Compressibility', 'Compressive Frequency': 'Peak Compressive Frequency', 'Transverse Frequency': 'Peak Transverse Frequency', 'Ellipticity': 'Ratio of Perpendicular Power', 'Takahashi Transverse Error/Resolution':'Takahashi Transverse Error/Measurement Resolution', 'Takahashi Compressive Error/Resolution':'Takahashi Compressive Error/Measurement Resolution', 'Heilig Transverse Error/Resolution':'Heilig Transverse Error/Measurement Resolution', 'Takahashi Transverse Difference': 'Takahashi Transverse Difference', 'Takahashi Compressive Difference': 'Takahashi Compressive Difference', 'Heilig Transverse Difference': 'Heilig Transverse Difference'}

    if 'obs_min' in choices:
        min_obs=choices['obs_min']
    else:
        min_obs=50

    hist_dict = {}
    
    #calculate basic histogram in all cases
    hist_count,x_edge, y_edge, _ = stats.binned_statistic_2d(df[x_col].to_numpy(), df[y_col].to_numpy(), df[w_dict['Transverse Power']].to_numpy(), statistic='count', bins=[x_bin_edges, y_bin_edges])

    #produce a copy of count distribution histogram for masking purposes
    hist_count_c = hist_count.copy()
    hist_count_c[hist_count_c < min_obs] = np.nan
    
    hist_count = hist_count.T
    hist_dict['count']=hist_count

    for arg in keys:
        arg_weight = w_dict[arg]
        arg_hist, _, _, _ = stats.binned_statistic_2d(df[x_col].to_numpy(), df[y_col].to_numpy(), df[arg_weight].to_numpy(), statistic=metric, bins=[x_bin_edges, y_bin_edges])
        arg_hist = np.where(np.isnan(hist_count_c), np.nan, arg_hist)
        arg_hist = arg_hist.T
        hist_dict[arg]=arg_hist
        
    return hist_dict

##################################################################

def compute_error_hists(df):
    """Compute 2D heatmaps of peak frequency, ignoring bins with <50 intervals"""
    x_col='GIPM X (OMNI mean)'
    y_col='GIPM Y (OMNI mean)'
    w_tak_trans ='Takahashi Transverse Error'
    w_tak_comp = 'Takahashi Compressive Error'
    w_heilig_trans = 'Heilig Transverse Error'

    x_bin_edges = range(0, 20)
    y_bin_edges = range(-20, 20)
    
    hist, _, _ = np.histogram2d(
        df[x_col].to_numpy(),
        df[y_col].to_numpy(),
        bins=[x_bin_edges, y_bin_edges]
    )
    hist = hist.T
    hist[hist == 0] = np.nan

    #produce a copy of count distribution histogram for masking purposes
    hist_count = hist.copy()
    hist_count[hist_count < 50] = np.nan

    #Takahashi transverse error histogram
    hist_tak_trans, _, _ = np.histogram2d(
        df[x_col].to_numpy(),
        df[y_col].to_numpy(),
        bins=[x_bin_edges, y_bin_edges],
        weights=df[w_tak_trans].to_numpy()
    )
    hist_tak_trans = hist_tak_trans.T
    hist_tak_trans[hist_tak_trans == 0] = np.nan
    #normalise to find averages
    hist_tak_trans = hist_tak_trans/hist_count

    #Takahashi compressive error power histogram
    hist_tak_comp, _, _ = np.histogram2d(
        df[x_col].to_numpy(),
        df[y_col].to_numpy(),
        bins=[x_bin_edges, y_bin_edges],
        weights=df[w_tak_comp].to_numpy()
    )
    hist_tak_comp = hist_tak_comp.T
    hist_tak_comp[hist_tak_comp == 0] = np.nan
    #normalise to find averages
    hist_tak_comp = hist_tak_comp/hist_count

    #Takahashi transverse error histogram
    hist_heilig_trans, _, _ = np.histogram2d(
        df[x_col].to_numpy(),
        df[y_col].to_numpy(),
        bins=[x_bin_edges, y_bin_edges],
        weights=df[w_heilig_trans].to_numpy()
    )
    hist_heilig_trans = hist_heilig_trans.T
    hist_heilig_trans[hist_heilig_trans == 0] = np.nan
    #normalise to find averages
    hist_heilig_trans = hist_heilig_trans/hist_count
    
    return hist, hist_tak_trans, hist_tak_comp, hist_heilig_trans
