from scipy import stats
import datetime as dt
import pandas as pd
import numpy as np

import matplotlib 
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.gridspec import GridSpec
from gipm_plotting_modules import model_calcs, draw_background, set_limits
from pyspedas.cotrans_tools.minvar import minvar

##################################################
#minimum variance analysis using pyspedas

#apply pyspedas code to cut-down window for whichever spacecraft we are comparing to

#######################################

def minvar_analysis(t_s: pd.DataFrame):
    """Input a time series and return the rotated frame basis vectors & eigenvalues"""
    bx = t_s["Bx_gse"].values
    by = t_s["By_gse"].values
    bz = t_s["Bz_gse"].values

    #vrot is the rotated time series for this interval
    #v gives new basis vectors, w/ max Vi=v[*,0], int Vi=v[*,1], etc
    #w gives eigenvalues
    vrot, v, w = minvar(np.column_stack([bx,by,bz]))
    print(f"Max/int/min eigenvalues: {w}")
    if w[1]/w[2] < 5:
        print("MVA propagation determination questionable")

    t_s['B_l']=vrot[:,0]
    t_s['B_m']=vrot[:,1]
    t_s['B_n']=vrot[:,2]
    return(v, w, t_s)

def minvar_rot(t_s: pd.DataFrame, v: np.Array):
    """Input a time series and mva basis vectors, return the time series rotated into the mva frame"""
    bx = t_s["Bx_gse"].values
    by = t_s["By_gse"].values
    bz = t_s["Bz_gse"].values
    
    gse_ts = np.column_stack([bx,by,bz])

    vrot = np.array([np.dot(row, v) for row in gse_ts])
    #vrot is the rotated time series for this interval

    t_s['B_l']=vrot[:,0]
    t_s['B_m']=vrot[:,1]
    t_s['B_n']=vrot[:,2]
    
    return(t_s)

##################################################

def max_corr_slide_1_min_gse(sc, start_time, comparison_df):
    """Slide in a +- 30s window to maximise gse correlation"""
    max_r_mag=0
    s_max_f_mag=-30
    max_r_x=0
    s_max_f_x=-30
    max_r_y=0
    s_max_f_y=-30
    max_r_z=0
    s_max_f_z=-30

    for s in np.arange(-30, 30, 0.2):
        comp_mag = comparison_df['B_mag'].to_numpy()
        comp_Bx = comparison_df['Bx_gse'].to_numpy()
        comp_By = comparison_df['By_gse'].to_numpy()
        comp_Bz = comparison_df['Bz_gse'].to_numpy()
    
        new_int_start=start_time + dt.timedelta(seconds=s)
        new_int_end=start_time+dt.timedelta(minutes=4)
        ts_df=df_dict[sc]['Time Series']
        ts_short[sc] = ts_df.loc[((ts_df.index >= new_int_start) & (ts_df.index < new_int_end))]
        mag = ts_short[sc]['B_mag'].to_numpy()
        bx = ts_short[sc]['Bx_gse'].to_numpy()
        by = ts_short[sc]['By_gse'].to_numpy()
        bz = ts_short[sc]['Bz_gse'].to_numpy()
        if mag.shape[0] > comp_mag.shape[0]:
            length = comp_mag.shape[0]
            mag = mag[:length]
            bx = bx[:length]
            by = by[:length]
            bz = bz[:length]
        if mag.shape[0] < comp_mag.shape[0]:
            length = mag.shape[0]
            comp_mag = comp_mag[:length]
            comp_Bx = comp_Bx[:length]
            comp_By = comp_By[:length]
            comp_Bz = comp_Bz[:length]
            
        r_mag,_ = stats.pearsonr(mag, comp_mag)
        r_x,_ = stats.pearsonr(bx, comp_Bx)
        r_y,_ = stats.pearsonr(by, comp_By)
        r_z,_ = stats.pearsonr(bz, comp_Bz)

        if abs(r_mag) > abs(max_r_mag):
            max_r_mag=r_mag
            s_max_f_mag=s
        if abs(r_x) > abs(max_r_x):
            max_r_x=r_x
            s_max_f_x=s
        if abs(r_y) > abs(max_r_y):
            max_r_y=r_y
            s_max_f_y=s
        if abs(r_z) > abs(max_r_z):
            max_r_z=r_z
            s_max_f_z=s    
    print(f"Max Bmag correlation: {max_r_mag} time: {start_time+dt.timedelta(seconds=s_max_f_mag)}.")
    print(f"Max Bx_gse correlation: {max_r_x} time: {start_time+dt.timedelta(seconds=s_max_f_x)}.")
    print(f"Max By_gse correlation: {max_r_y} time: {start_time+dt.timedelta(seconds=s_max_f_y)}.")
    print(f"Max Bz_gse correlation: {max_r_z} time: {start_time+dt.timedelta(seconds=s_max_f_z)}.")

##################################################

def max_corr_slide_1_min_MVA(ts_df, start_time, comparison_df):
    """Slide in a +- 30s window to maximise mva correlation and then plot the two with the highest transverse correlation"""
    #lmn coordinate system (max/int/min var directions)
    max_r_mag=0
    s_max_f_mag=-30
    max_r_l=0
    s_max_f_l=-30
    max_r_m=0
    s_max_f_m=-30
    max_r_n=0
    s_max_f_n=-30

    for s in np.arange(-30, 30, 0.2):
        comp_mag = comparison_df['B_mag'].to_numpy()
        comp_Bl = comparison_df['B_l'].to_numpy()
        comp_Bm = comparison_df['B_m'].to_numpy()
        comp_Bn = comparison_df['B_n'].to_numpy()
    
        new_int_start=start_time + dt.timedelta(seconds=s)
        new_int_end=start_time+dt.timedelta(minutes=4)
        ts_short = ts_df.loc[((ts_df.index >= new_int_start) & (ts_df.index < new_int_end))]
        mag = ts_short['B_mag'].to_numpy()
        bl = ts_short['B_l'].to_numpy()
        bm = ts_short['B_m'].to_numpy()
        bn = ts_short['B_n'].to_numpy()
        if mag.shape[0] > comp_mag.shape[0]:
            length = comp_mag.shape[0]
            mag = mag[:length]
            bl = bl[:length]
            bm = bm[:length]
            bn = bn[:length]
        if mag.shape[0] < comp_mag.shape[0]:
            length = mag.shape[0]
            comp_mag = comp_mag[:length]
            comp_Bl = comp_Bl[:length]
            comp_Bm = comp_Bm[:length]
            comp_Bn = comp_Bn[:length]
            
        r_mag,_ = stats.pearsonr(mag, comp_mag)
        r_l,_ = stats.pearsonr(bl, comp_Bl)
        r_m,_ = stats.pearsonr(bm, comp_Bm)
        r_n,_ = stats.pearsonr(bn, comp_Bn)

        if abs(r_mag) > abs(max_r_mag):
            max_r_mag=r_mag
            s_max_f_mag=s
        if abs(r_l) > abs(max_r_l):
            max_r_l=r_l
            s_max_f_l=s
        if abs(r_m) > abs(max_r_m):
            max_r_m=r_m
            s_max_f_m=s
        if abs(r_n) > abs(max_r_n):
            max_r_n=r_n
            s_max_f_n=s    

    mr_mag_t=start_time+dt.timedelta(seconds=s_max_f_mag)
    mr_l_t=start_time+dt.timedelta(seconds=s_max_f_l)
    mr_m_t=start_time+dt.timedelta(seconds=s_max_f_m)
    mr_n_t=start_time+dt.timedelta(seconds=s_max_f_n)
    print(f"Max Bmag correlation: {max_r_mag} time: {mr_mag_t}.")
    print(f"Max B_l correlation: {max_r_l} time: {mr_l_t}.")
    print(f"Max B_m correlation: {max_r_m} time: {mr_m_t}.")
    print(f"Max B_n correlation: {max_r_n} time: {mr_n_t}.")

    #return time of max bmag/l/m,n correlation
    return(mr_mag_t, mr_l_t, mr_m_t, mr_n_t)

##################################################

def ts_overlay(inp_ts_1: pd.DataFrame, inp_ts_2: pd.DataFrame, fname: str, scs: list[str], props=['B_mag']) -> None: 
    """Overlay time series for input time series"""
    fig = plt.figure(figsize=(8, (len(props)//2)*2))

    if len(props)==1:
        ax=fig.add_axes((0, 0, 1, 1))
        ax_list=[ax]
    else:
        ax_list=[]
        gs = GridSpec(len(props)//2, 2, figure=fig)
        for n in range(len(props)):
            ax = fig.add_subplot(gs[n//2,n%2])
            ax_list.append(ax)
    
    locator = mdates.AutoDateLocator(minticks=3, maxticks=5)
    formatter = mdates.ConciseDateFormatter(locator)

    for l,p in enumerate(props):
        ax=ax_list[l]
        if inp_ts_1[p].min() < inp_ts_2[p].min():
            Cmin=(inp_ts_1[p].min() - 1)
        else:
            Cmin=(inp_ts_2[p].min() - 1)
        if inp_ts_1[p].max() > inp_ts_2[p].max():
            Cmax=(inp_ts_1[p].max()+ 1)
        else:
            Cmax=(inp_ts_2[p].max() + 1)        
        ax.plot(inp_ts_1.index, inp_ts_1[p], color = 'k', label=scs[0])
        ax.plot(inp_ts_1.index, inp_ts_2[p], color = 'indigo', label=scs[1])
        ax.set_ylim(Cmin, Cmax)
        ax.set_xlim(inp_ts_1.index[0], inp_ts_1.index[-1])
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.xaxis.set_tick_params(labelsize=12)
        ax.yaxis.set_tick_params(labelsize=12)
        ax.set_ylabel(p)
        if p==props[-1]:
            ax.legend(loc='lower right')
    
    path = f"/Users/roseatkinson/Documents/New_Figs/Correlation_Study/{fname}.png"
    plt.savefig(path, bbox_inches='tight')

##################################################

def time_series_comp(df_dict:dict, sc_list:list[str],fname:str) -> None:
    """Take time series dict and spacecraft reference to plot two ts side-by-side"""
    if len(sc_list) > 2:
        print("Spacecraft list too long, only handling first two")
        sc_list=sc_list[:2]
    
    fig = plt.figure(figsize=(8, 3))
        
    locator = mdates.AutoDateLocator(minticks=3, maxticks=5)
    formatter = mdates.ConciseDateFormatter(locator)
        
    # Define a gridspec layout
    gs = GridSpec(1, 2, figure=fig)
    
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    
    for ax, sc in zip([ax1,ax2], sc_list):
        df=df_dict[sc]
        ax.plot(df.index, df['Bx_gse'], color = 'red', label='Bx')
        ax.plot(df.index, df['By_gse'], color = 'green', label='By')
        ax.plot(df.index, df['Bz_gse'], color = 'blue', label='Bz')
        ax.plot(df.index, df['B_mag'], color = 'k', label='Bmag')
        ax.set_ylim(-10, 10)
        ax.set_xlim(df.index[0], df.index[-1])
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.xaxis.set_tick_params(labelsize=12)
        ax.yaxis.set_tick_params(labelsize=12)
        ax.set_ylabel(r'B, nT')
        ax.set_title(sc)
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='lightgray')
        ax.xaxis.grid(color='lightgray')
        if sc==sc_list[1]:
            ax.legend(loc='upper right')
    
    path = f"/Users/roseatkinson/Documents/New_Figs/Correlation_Study/{fname}.png"
            
    plt.savefig(path, bbox_inches='tight')

##################################################

def loc_plot(locs_x:list[float], locs_y:list[float], fname:str) -> None: 
    
    X_shue, R_shue, Xgipm, Ygipm, Zgipm, f, fitting_coeffs = model_calcs()

    """Plots two locations in X,Y GIPM"""
    fig = plt.figure(figsize=(4, 4))
    ax=fig.add_axes((0, 0, 1, 1))
    
    draw_background(ax, Xgipm[:, :, 0], Ygipm[:, :, 0], f[:, :, 0],
                        X_shue, R_shue)
    
    ax.scatter(locs_x, locs_y, marker="o", c=['k','indigo'], s=60)
    
    set_limits(ax)
    
    ax.set_title(r'Case Study Locations')
    ax.set_xlabel("$X_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")
    ax.set_ylabel("$Y_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")
    
    path = f"/Users/roseatkinson/Documents/New_Figs/Correlation_Study/{fname}.png"
    plt.savefig(path, bbox_inches='tight')

