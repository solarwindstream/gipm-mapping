def compute_hists2d(df):
    """Compute 2D heatmaps of power and compressibilty, ignoring bins with <50 intervals"""
    x_col='GIPM X (OMNI mean)'
    y_col='GIPM Y (OMNI mean)'
    w_compressive='ULF Band Normalised Compressive Power'
    w_transverse='ULF Band Normalised Transverse Power'
    w_compressible='Compressibility'
    
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

    #compressibility histogram
    hist_compressibility, _, _ = np.histogram2d(
        df[x_col].to_numpy(),
        df[y_col].to_numpy(),
        bins=[x_bin_edges, y_bin_edges],
        weights=df[w_compressible].to_numpy()
    )
    hist_compressibility = hist_compressibility.T
    hist_compressibility[hist_compressibility == 0] = np.nan
    #normalise to find averages
    hist_compressibility = hist_compressibility/hist_count
    
    return hist, hist_comp, hist_trans, hist_compressibility


def compute_freq_ellip_hists(df):
    """Compute 2D heatmaps of peak frequency and ellipticity, ignoring bins with <50 intervals"""
    x_col='GIPM X (OMNI mean)'
    y_col='GIPM Y (OMNI mean)'
    w_compressive='Peak Compressive Frequency'
    w_transverse='Peak Transverse Frequency'
    w_ellipticity='Ratio of Perpendicular Power'
    
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