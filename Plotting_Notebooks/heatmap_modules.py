import numpy as np 
import pandas as pd
from scipy import stats

##################################################################

def compute_hists2d(df):
    """Compute 2D heatmaps of transverse and compressive power and compressibilty, ignoring bins with <50 intervals"""
    x_col='GIPM X (OMNI mean)'
    y_col='GIPM Y (OMNI mean)'
    w_transverse='ULF Band Normalised Transverse Power'
    w_compressive='ULF Band Normalised Compressive Power'
    w_compressible='Compressibility'

    hist_count,x_edge, y_edge, _ = stats.binned_statistic_2d(df[x_col].to_numpy(), df[y_col].to_numpy(), df[w_transverse].to_numpy(), statistic='count', bins=[x_bin_edges, y_bin_edges])
    
    #produce a copy of count distribution histogram for masking purposes
    hist_count_c = hist.copy()
    hist_count_c[hist_count < 50] = np.nan

    #transverse
    hist_trans,x_edge, y_edge, _ = stats.binned_statistic_2d(df[x_col].to_numpy(), df[y_col].to_numpy(), df[w_transverse].to_numpy(), statistic='mean', bins=[x_bin_edges, y_bin_edges])
    hist_trans = np.where(np.isnan(hist_count_c), np.nan, hist_trans)
    hist_trans = hist_trans.T

    #transverse std
    hist_trans_std,x_edge, y_edge, _ = stats.binned_statistic_2d(df[x_col].to_numpy(), df[y_col].to_numpy(), df[w_transverse].to_numpy(), statistic='std', bins=[x_bin_edges, y_bin_edges])
    hist_trans_std = np.where(np.isnan(hist_count_c), np.nan, hist_trans_std)
    hist_trans_std = hist_trans_std.T

    #compressive
    hist_comp,x_edge, y_edge, _ = stats.binned_statistic_2d(df[x_col].to_numpy(), df[y_col].to_numpy(), df[w_compressive].to_numpy(), statistic='mean', bins=[x_bin_edges, y_bin_edges])
    hist_comp = np.where(np.isnan(hist_count_c), np.nan, hist_comp)
    hist_comp = hist_comp.T

    #compressive std
    hist_comp_std,x_edge, y_edge, _ = stats.binned_statistic_2d(df[x_col].to_numpy(), df[y_col].to_numpy(), df[w_compressive].to_numpy(), statistic='std', bins=[x_bin_edges, y_bin_edges])
    hist_comp_std = np.where(np.isnan(hist_count_c), np.nan, hist_comp_std)
    hist_comp_std = hist_comp_std.T
    
    #compressibility
    hist_compressibility,x_edge, y_edge, _ = stats.binned_statistic_2d(df[x_col].to_numpy(), df[y_col].to_numpy(), df[w_compressible].to_numpy(), statistic='mean', bins=[x_bin_edges, y_bin_edges])
    hist_compressibility = np.where(np.isnan(hist_count_c), np.nan, hist_compressibility)
    hist_compressibility = hist_compressibility.T

    hist_count = hist_count.T
    print("results order is hist/transverse mean/trans std/comp mean/comp std/compressibility mean now! update me!")
    return hist_count, hist_trans, hist_trans_std, hist_comp, hist_comp_std, hist_compressibility

##################################################################

def compute_hists2d_low_data(df):
    """Compute 2D heatmaps of power, ignoring bins with <30 intervals"""
    x_col='GIPM X (OMNI mean)'
    y_col='GIPM Y (OMNI mean)'
    w_compressive='ULF Band Normalised Compressive Power'
    w_transverse='ULF Band Normalised Transverse Power'

    x_bin_edges = range(20)
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
    hist_count[hist_count < 20] = np.nan

    #compressive power histogram
    hist_comp, _, _ = np.histogram2d(
        df[x_col].to_numpy(),
        df[y_col].to_numpy(),
        bins=[x_bin_edges, y_bin_edges],
        weights=df[w_compressive].to_numpy()
    )
    hist_comp = hist_comp.T
    hist_comp[hist_comp == 0] = np.nan
    #normalise to find averages
    hist_comp = hist_comp/hist_count

    #transverse power histogram
    hist_trans, _, _ = np.histogram2d(
        df[x_col].to_numpy(),
        df[y_col].to_numpy(),
        bins=[x_bin_edges, y_bin_edges],
        weights=df[w_transverse].to_numpy()
    )
    hist_trans = hist_trans.T
    hist_trans[hist_trans == 0] = np.nan
    #normalise to find averages
    hist_trans = hist_trans/hist_count

    return hist, hist_comp, hist_trans

##################################################################

def compute_freq_ellip_hists(df):
    """Compute 2D heatmaps of peak frequency and ellipticity, ignoring bins with <50 intervals"""
    x_col='GIPM X (OMNI mean)'
    y_col='GIPM Y (OMNI mean)'
    w_compressive='Peak Compressive Frequency'
    w_transverse='Peak Transverse Frequency'
    w_ellipticity='Ratio of Perpendicular Power'

    x_bin_edges = range(20)
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

    #compressive power histogram
    hist_comp, _, _ = np.histogram2d(
        df[x_col].to_numpy(),
        df[y_col].to_numpy(),
        bins=[x_bin_edges, y_bin_edges],
        weights=df[w_compressive].to_numpy()
    )
    hist_comp = hist_comp.T
    hist_comp[hist_comp == 0] = np.nan
    #normalise to find averages
    hist_comp = hist_comp/hist_count

    #transverse power histogram
    hist_trans, _, _ = np.histogram2d(
        df[x_col].to_numpy(),
        df[y_col].to_numpy(),
        bins=[x_bin_edges, y_bin_edges],
        weights=df[w_transverse].to_numpy()
    )
    hist_trans = hist_trans.T
    hist_trans[hist_trans == 0] = np.nan
    #normalise to find averages
    hist_trans = hist_trans/hist_count

    #ellipticity histogram
    hist_ellipticity, _, _ = np.histogram2d(
        df[x_col].to_numpy(),
        df[y_col].to_numpy(),
        bins=[x_bin_edges, y_bin_edges],
        weights=df[w_ellipticity].to_numpy()
    )
    hist_ellipticity = hist_ellipticity.T
    hist_ellipticity[hist_ellipticity == 0] = np.nan
    #normalise to find averages
    hist_ellipticity = hist_ellipticity/hist_count
    
    return hist, hist_comp, hist_trans, hist_ellipticity

##################################################################

def compute_normalised_freq_hists(df):
    """Compute 2D heatmaps of peak *normalised* frequency, ignoring bins with <50 intervals"""
    x_col='GIPM X (OMNI mean)'
    y_col='GIPM Y (OMNI mean)'
    w_compressive='Normalised Compressive Frequency'
    w_transverse='Normalised Transverse Frequency'

    x_bin_edges = range(20)
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

    #compressive power histogram
    hist_comp, _, _ = np.histogram2d(
        df[x_col].to_numpy(),
        df[y_col].to_numpy(),
        bins=[x_bin_edges, y_bin_edges],
        weights=df[w_compressive].to_numpy()
    )
    hist_comp = hist_comp.T
    hist_comp[hist_comp == 0] = np.nan
    #normalise to find averages
    hist_comp = hist_comp/hist_count

    #transverse power histogram
    hist_trans, _, _ = np.histogram2d(
        df[x_col].to_numpy(),
        df[y_col].to_numpy(),
        bins=[x_bin_edges, y_bin_edges],
        weights=df[w_transverse].to_numpy()
    )
    hist_trans = hist_trans.T
    hist_trans[hist_trans == 0] = np.nan
    #normalise to find averages
    hist_trans = hist_trans/hist_count
    
    return hist, hist_comp, hist_trans

##################################################################
    
def compute_freq_hists_low_data(df):
    """Compute 2D heatmaps of peak frequency, ignoring bins with <20 intervals"""
    x_col='GIPM X (OMNI mean)'
    y_col='GIPM Y (OMNI mean)'
    w_compressive='Peak Compressive Frequency'
    w_transverse='Peak Transverse Frequency'

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

    #compressive power histogram
    hist_comp, _, _ = np.histogram2d(
        df[x_col].to_numpy(),
        df[y_col].to_numpy(),
        bins=[x_bin_edges, y_bin_edges],
        weights=df[w_compressive].to_numpy()
    )
    hist_comp = hist_comp.T
    hist_comp[hist_comp == 0] = np.nan
    #normalise to find averages
    hist_comp = hist_comp/hist_count

    #transverse power histogram
    hist_trans, _, _ = np.histogram2d(
        df[x_col].to_numpy(),
        df[y_col].to_numpy(),
        bins=[x_bin_edges, y_bin_edges],
        weights=df[w_transverse].to_numpy()
    )
    hist_trans = hist_trans.T
    hist_trans[hist_trans == 0] = np.nan
    #normalise to find averages
    hist_trans = hist_trans/hist_count
    
    return hist, hist_comp, hist_trans


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

######################

def compute_res_error_hists(df):
    """Compute 2D heatmaps of peak frequency, ignoring bins with <50 intervals"""
    x_col='GIPM X (OMNI mean)'
    y_col='GIPM Y (OMNI mean)'
    w_tak_trans ='Takahashi Transverse Error/Measurement Resolution'
    w_tak_comp ='Takahashi Compressive Error/Measurement Resolution'
    w_heilig_trans ='Heilig Transverse Error/Measurement Resolution'
    
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

    #Takahashi compressive error histogram
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

    #Heilig transverse error histogram
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
