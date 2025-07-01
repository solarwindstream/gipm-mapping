#Plot Hann Windowed Time Series

import pandas as pd
import numpy as np
from numpy.fft import fft, ifft, rfft
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.ticker as ticker
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.dates as mdates

pd.options.mode.chained_assignment = None

def Hann_TS(cluster_ULF_csv, str_centre):
    
    #dl_freq = dl_freq * (10**(-3))
    #tak_freq = tak_freq * (10**(-3))
    
    time_cent = pd.to_datetime(str_centre)
    
    ULF_df = pd.read_csv(cluster_ULF_csv)
    ULF_df['datetime'] = pd.to_datetime(ULF_df['datetime'])
    ULF_df = ULF_df.set_index('datetime')

    locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
    formatter = mdates.ConciseDateFormatter(locator)

    #mask to just relevant times
    
    time_start_4mins = time_cent - dt.timedelta(minutes=2)
    time_end_4mins = time_cent + dt.timedelta(minutes=2)

    ULF_df_4mins = ULF_df.loc[((ULF_df.index >= time_start_4mins) & (ULF_df.index < time_end_4mins))]
    
    # sampling rate
    sr = 22
    sample_rate = 1/22
    ecf = np.sqrt(8/3)
    
   
    #######FOUR MINUTES
    
    #now, find average Cluster magnetic field direction during this time

    #Bx_gse, By_gse, Bz_gse, B_mag

    #normalised columns for each, then average

    ULF_df_4mins['Norm_Bx'] = ULF_df_4mins['Bx_gse'].div(ULF_df_4mins['B_mag'])
    ULF_df_4mins['Norm_By'] = ULF_df_4mins['By_gse'].div(ULF_df_4mins['B_mag'])
    ULF_df_4mins['Norm_Bz'] = ULF_df_4mins['Bz_gse'].div(ULF_df_4mins['B_mag'])

    mean_x = ULF_df_4mins['Norm_Bx'].mean()
    mean_y = ULF_df_4mins['Norm_By'].mean()
    mean_z = ULF_df_4mins['Norm_Bz'].mean()

    norm = (mean_x**2 + mean_y**2 + mean_z**2)**0.5

    mean_x = mean_x/norm
    mean_y = mean_y/norm
    mean_z = mean_z/norm

    #make a new basis w/ Mean B, vector perp to Mean B and X, and third perpendicular vector.
    ##NB NORMALISE

    B_mean = np.array([mean_x, mean_y, mean_z])
    X_dir = np.array([1,0,0])

    Norm_1 = np.cross(B_mean, X_dir)
    Norm_1 = Norm_1/(np.linalg.norm(Norm_1))
    Norm_2 = np.cross(B_mean, Norm_1)
    Norm_2 = Norm_2/(np.linalg.norm(Norm_2))
    #find B in direction of B perp 
    #first make array for each row in table, then get dot products for each 
    
    print('B mean 4mins:',   B_mean)
    print('Norm 1 4mins:',  Norm_1)
    print('Norm 2 4mins:',  Norm_2)
    
    list_of_B = []

    for i,j,k in zip(ULF_df_4mins['Bx_gse'], ULF_df_4mins['By_gse'], ULF_df_4mins['Bz_gse']):
        B_vec = np.array([i,j,k])
        list_of_B.append(B_vec)

    ULF_df_4mins['B Vector']= list_of_B

    list_norm_1_vals = []

    for i in ULF_df_4mins['B Vector']:
        projection_norm_1 = np.dot(i, Norm_1)
        list_norm_1_vals.append(projection_norm_1)

    ULF_df_4mins['B Perp 1'] = list_norm_1_vals

    list_norm_2_vals = []

    for i in ULF_df_4mins['B Vector']:
        projection_norm_2 = np.dot(i, Norm_2)
        list_norm_2_vals.append(projection_norm_2)
        
    list_para_vals = []

    for i in ULF_df_4mins['B Vector']:
        projection_para = np.dot(i, B_mean)
        list_para_vals.append(projection_para)

    ULF_df_4mins['B Perp 2'] = list_norm_2_vals
    ULF_df_4mins['B Para'] = list_para_vals
    
    para_mean_4 = ULF_df_4mins['B Para'].mean()
    ULF_df_4mins['B Para'] = ULF_df_4mins['B Para'] - para_mean_4

    #now calculate power spectrum for each perp component & sum.

    #FFT 2 mins norm 1
    x_4mins_perp_1 = ULF_df_4mins['B Perp 1'].to_numpy()
    four_min_hann = np.hanning(len(x_4mins_perp_1))
    x_4mins_perp_1 = x_4mins_perp_1*four_min_hann

    X_4mins_perp_1 = rfft(x_4mins_perp_1,norm='ortho')
    N_4mins_perp_1 = len(X_4mins_perp_1)
    n_4mins_perp_1 = np.arange(N_4mins_perp_1)
    T_4mins_perp_1 = N_4mins_perp_1/sr
    freq_4mins_perp_1 = n_4mins_perp_1/T_4mins_perp_1
    power_4mins_perp_1 = (np.abs(X_4mins_perp_1)**2)*ecf*sample_rate

    #FFT 2 mins norm 2
    x_4mins_perp_2 = ULF_df_4mins['B Perp 2'].to_numpy()
    x_4mins_perp_2 = x_4mins_perp_2*four_min_hann

    X_4mins_perp_2 = rfft(x_4mins_perp_2,norm='ortho')
    N_4mins_perp_2 = len(X_4mins_perp_2)
    n_4mins_perp_2 = np.arange(N_4mins_perp_2)
    T_4mins_perp_2 = N_4mins_perp_2/sr
    freq_4mins_perp_2 = n_4mins_perp_2/T_4mins_perp_2
    power_4mins_perp_2 = (np.abs(X_4mins_perp_2)**2)*ecf*sample_rate


    #add up (since power, do not need to sqrt). wait but freqs might not be same. 

    power_4mins_perp = power_4mins_perp_1 + power_4mins_perp_2

    #power_tot = power_2_1 + power_2_2
    
    #FFT parallel
    x_4mins_para = ULF_df_4mins['B Para'].to_numpy()
    x_4mins_para = x_4mins_para*four_min_hann

    X_4mins_para = rfft(x_4mins_para,norm='ortho')
    N_4mins_para = len(X_4mins_para)
    n_4mins_para = np.arange(N_4mins_para)
    T_4mins_para = N_4mins_para/sr
    freq_4mins_para = n_4mins_para/T_4mins_para
    power_4mins_para = (np.abs(X_4mins_para)**2)*ecf*sample_rate
    
    #show windowed vs unwindowed components
    
    fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 8), tight_layout=True)
    
    #fig.suptitle('Interval Centred on ' + str_centre)
    ax1.set_title('Raw Data, Mean Removed')
    ax1.plot(ULF_df_4mins.index, ULF_df_4mins['B Para'], color = 'g', label='$B_{\parallel}$')
    ax1.plot(ULF_df_4mins.index, ULF_df_4mins['B Perp 1'], color = 'b', label='$B_{\perp1}$')
    ax1.plot(ULF_df_4mins.index, ULF_df_4mins['B Perp 2'], color = 'r', label='$B_{\perp2}$')
    ax1.plot(ULF_df_4mins.index, four_min_hann, color = 'k', label='Hann Window')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Magnetic Field, nT')
    ax1.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    ax1.set_axisbelow(True)
    ax1.xaxis.set_major_formatter(formatter)
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.grid(color='lightgray')
    ax1.xaxis.grid(color='lightgray')
    
    
    ax2.set_title('Hann Window Applied')
    ax2.plot(ULF_df_4mins.index, x_4mins_para, color = 'g', label='$B_{\parallel}$')
    ax2.plot(ULF_df_4mins.index, x_4mins_perp_1, color = 'b', label='$B_{\perp1}$')
    ax2.plot(ULF_df_4mins.index, x_4mins_perp_2, color = 'r', label='$B_{\perp2}$')
    ax2.plot(ULF_df_4mins.index, four_min_hann, color = 'k', label='Hann Window')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Magnetic Field, nT')
    ax2.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    ax2.xaxis.set_major_formatter(formatter)
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.set_axisbelow(True)
    ax2.yaxis.grid(color='lightgray')
    ax2.xaxis.grid(color='lightgray')
    
    fig.show()