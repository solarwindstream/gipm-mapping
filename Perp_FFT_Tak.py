#Fourier Transform B Perp edition

import pandas as pd
import numpy as np
from numpy.fft import fft, ifft, rfft
import matplotlib.pyplot as plt
import datetime as dt

def FFT_perp_Tak(cluster_ULF_csv, str_centre, dl_freq, tak_freq):
    
    dl_freq = dl_freq * (10**(-3))
    tak_freq = tak_freq * (10**(-3))
    
    time_cent = pd.to_datetime(str_centre)
    
    ULF_df = pd.read_csv(cluster_ULF_csv)
    ULF_df['datetime'] = pd.to_datetime(ULF_df['datetime'])
    ULF_df = ULF_df.set_index('datetime')

    #mask to just relevant times
    #time_start_1min = time_cent - dt.timedelta(seconds=30)
    #time_end_1min = time_cent + dt.timedelta(seconds=30)

    time_start_2mins = time_cent - dt.timedelta(minutes=1)
    time_end_2mins = time_cent + dt.timedelta(minutes=1)

    #time_start_2mins = pd.to_datetime('2001-04-23 05:15:00')
    #time_end_2mins = pd.to_datetime('2001-04-23 05:19:00')

    #time_start_10mins = pd.to_datetime('2001-04-23 05:12:00')
    #time_end_10mins = pd.to_datetime('2001-04-23 05:22:00')

    #ULF_df_1min = ULF_df.loc[((ULF_df.index >= time_start_1min) & (ULF_df.index < time_end_1min))]
    ULF_df_2mins = ULF_df.loc[((ULF_df.index >= time_start_2mins) & (ULF_df.index < time_end_2mins))]
    #ULF_df_4mins = ULF_df.loc[((ULF_df.index >= time_start_4mins) & (ULF_df.index < time_end_4mins))]
    #ULF_df_10mins = ULF_df.loc[((ULF_df.index >= time_start_10mins) & (ULF_df.index < time_end_10mins))]

    #now, find average Cluster magnetic field direction during this time

    #Bx_gse, By_gse, Bz_gse, B_mag

    #normalised columns for each, then average

    ULF_df_2mins['Norm_Bx'] = ULF_df_2mins['Bx_gse']/ULF_df_2mins['B_mag']
    ULF_df_2mins['Norm_By'] = ULF_df_2mins['By_gse']/ULF_df_2mins['B_mag']
    ULF_df_2mins['Norm_Bz'] = ULF_df_2mins['Bz_gse']/ULF_df_2mins['B_mag']

    mean_x = ULF_df_2mins['Norm_Bx'].mean()
    mean_y = ULF_df_2mins['Norm_By'].mean()
    mean_z = ULF_df_2mins['Norm_Bz'].mean()

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
    Norm_2 = np.cross(X_dir, Norm_1)
    Norm_2 = Norm_2/(np.linalg.norm(Norm_2))
    #find B in direction of B perp 
    #first make array for each row in table, then get dot products for each 

    list_of_B = []

    for i,j,k in zip(ULF_df_2mins['Bx_gse'], ULF_df_2mins['By_gse'], ULF_df_2mins['Bz_gse']):
        B_vec = np.array([i,j,k])
        list_of_B.append(B_vec)

    ULF_df_2mins['B Vector']= list_of_B

    list_norm_1_vals = []

    for i in ULF_df_2mins['B Vector']:
        projection_norm_1 = np.dot(i, Norm_1)
        list_norm_1_vals.append(projection_norm_1)

    ULF_df_2mins['B Perp 1'] = list_norm_1_vals

    list_norm_2_vals = []

    for i in ULF_df_2mins['B Vector']:
        projection_norm_2 = np.dot(i, Norm_2)
        list_norm_2_vals.append(projection_norm_2)
        
    list_para_vals = []

    for i in ULF_df_2mins['B Vector']:
        projection_para = np.dot(i, B_mean)
        list_para_vals.append(projection_para)

    ULF_df_2mins['B Perp 2'] = list_norm_2_vals
    ULF_df_2mins['B Para'] = list_para_vals

    #now calculate power spectrum for each perp component & sum.

    #FFT 2 mins norm 1
    x_2_1 = ULF_df_2mins['B Perp 1'].to_numpy()
    # sampling rate
    sr = 22

    X_2_1 = rfft(x_2_1)
    N_2_1 = len(X_2_1)
    n_2_1 = np.arange(N_2_1)
    T_2_1 = N_2_1/sr
    freq_2_1 = n_2_1/T_2_1
    power_2_1 = np.abs(X_2_1)**2

    #FFT 2 mins norm 2
    x_2_2 = ULF_df_2mins['B Perp 2'].to_numpy()
    # sampling rate
    sr = 22

    X_2_2 = rfft(x_2_2)
    N_2_2 = len(X_2_2)
    n_2_2 = np.arange(N_2_2)
    T_2_2 = N_2_2/sr
    freq_2_2 = n_2_2/T_2_2
    power_2_2 = np.abs(X_2_2)**2


    #add up (since power, do not need to sqrt). wait but freqs might not be same. 

    power_2_t = power_2_1 + power_2_2

    #power_tot = power_2_1 + power_2_2
    
    #FFT parallel
    x_2_3 = ULF_df_2mins['B Para'].to_numpy()
    # sampling rate
    sr = 22

    X_2_3 = rfft(x_2_3)
    N_2_3 = len(X_2_3)
    n_2_3 = np.arange(N_2_3)
    T_2_3 = N_2_3/sr
    freq_2_3 = n_2_3/T_2_3
    power_2_3 = np.abs(X_2_3)**2
    
    
    #################

    #show spectra in two perp directions.
    fig, (ax1) = plt.subplots(1)
    plt.suptitle(str(time_start_2mins)+'-'+str(time_end_2mins))
    ax1.plot(freq_2_3, power_2_3, color="red", label='Parallel Power')
    ax1.plot(freq_2_1, power_2_t, color="black", label='Perpendicular Power')
    #plt.xlabel('Period (s)')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel('Power')
    ax1.set_xlabel('Frequency, Hz')
    ax1.set_xlim(0.01, 5)
    ax1.set_ylim(0.1, 100_000_000)
    ax1.vlines(x=dl_freq, ymin = 0.1, ymax = 100_000_000, linestyles='dashed', label='De Lauretis Freq')
    ax1.vlines(x=tak_freq, ymin = 0.1, ymax = 100_000_000, linestyles='dotted', label='Takahashi Freq')
    ax1.legend()

    plt.tight_layout()
    plt.show()