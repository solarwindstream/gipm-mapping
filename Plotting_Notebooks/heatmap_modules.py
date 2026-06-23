import numpy as np 
import pandas as pd
from scipy import stats

##################################################################

def compute_hists2d(df, metric, *keys, **choices):
    """Compute 2D heatmaps of transverse/compressive power/frequency/compressibilty as specified, ignoring bins with obs_min no. of intervals (50 default)"""

    x_bin_edges = range(20)
    y_bin_edges = range(-20, 20)
    
    x_col='GIPM X (OMNI mean)'
    y_col='GIPM Y (OMNI mean)'

    w_dict = {'Transverse Power':'ULF Band Normalised Transverse Power', 'Compressive Power':'ULF Band Normalised Compressive Power', 'Compressibility':'Compressibility', 'Compressive Frequency': 'Peak Compressive Frequency', 'Transverse Frequency': 'Peak Transverse Frequency', 'Ellipticity': 'Ratio of Perpendicular Power', 'Takahashi Transverse Error/Resolution':'Takahashi Transverse Error/Measurement Resolution', 'Takahashi Compressive Error/Resolution':'Takahashi Compressive Error/Measurement Resolution', 'Heilig Transverse Error/Resolution':'Heilig Transverse Error/Measurement Resolution', 'Takahashi Transverse Difference': 'Takahashi Transverse Difference', 'Takahashi Compressive Difference': 'Takahashi Compressive Difference', 'Heilig Transverse Difference': 'Heilig Transverse Difference', 'Heilig Compressive Difference': 'Heilig Compressive Difference'}

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
