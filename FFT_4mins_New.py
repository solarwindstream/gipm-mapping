#Fourier Transform B Perp edition

import pandas as pd
import numpy as np
from numpy.fft import fft, ifft, rfft
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.ticker as ticker

pd.options.mode.chained_assignment = None

def FFT_perp_4(cluster_ULF_csv, str_centre):
    
    #dl_freq = dl_freq * (10**(-3))
    #tak_freq = tak_freq * (10**(-3))
    
    time_cent = pd.to_datetime(str_centre)
    
    ULF_df = pd.read_csv(cluster_ULF_csv)
    ULF_df['datetime'] = pd.to_datetime(ULF_df['datetime'])
    ULF_df = ULF_df.set_index('datetime')

    #mask to just relevant times

    time_start_2mins = time_cent - dt.timedelta(minutes=1)
    time_end_2mins = time_cent + dt.timedelta(minutes=1)
    time_start_4mins = time_cent - dt.timedelta(minutes=2)
    time_end_4mins = time_cent + dt.timedelta(minutes=2)
    
    ULF_df_2mins = ULF_df.loc[((ULF_df.index >= time_start_2mins) & (ULF_df.index < time_end_2mins))]
    ULF_df_4mins = ULF_df.loc[((ULF_df.index >= time_start_4mins) & (ULF_df.index < time_end_4mins))]
    
    # sampling rate
    sr = 22
    sample_rate = 1/22
    
    int_lower_lim = 7*(10**(-3))
    int_upper_lim = 100*(10**(-3))
    
   
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

    X_4mins_perp_1 = rfft(x_4mins_perp_1,norm='ortho')
    N_4mins_perp_1 = len(X_4mins_perp_1)
    n_4mins_perp_1 = np.arange(N_4mins_perp_1)
    T_4mins_perp_1 = N_4mins_perp_1/sr
    freq_4mins_perp_1 = n_4mins_perp_1/T_4mins_perp_1
    power_4mins_perp_1 = (np.abs(X_4mins_perp_1)**2)*sample_rate

    #FFT 2 mins norm 2
    x_4mins_perp_2 = ULF_df_4mins['B Perp 2'].to_numpy()

    X_4mins_perp_2 = rfft(x_4mins_perp_2,norm='ortho')
    N_4mins_perp_2 = len(X_4mins_perp_2)
    n_4mins_perp_2 = np.arange(N_4mins_perp_2)
    T_4mins_perp_2 = N_4mins_perp_2/sr
    freq_4mins_perp_2 = n_4mins_perp_2/T_4mins_perp_2
    power_4mins_perp_2 = (np.abs(X_4mins_perp_2)**2)*sample_rate


    #add up (since power, do not need to sqrt)

    power_4mins_perp = power_4mins_perp_1 + power_4mins_perp_2

    #power_tot = power_2_1 + power_2_2
    
    #FFT parallel
    x_4mins_para = ULF_df_4mins['B Para'].to_numpy()

    X_4mins_para = rfft(x_4mins_para,norm='ortho')
    N_4mins_para = len(X_4mins_para)
    n_4mins_para = np.arange(N_4mins_para)
    T_4mins_para = N_4mins_para/sr
    freq_4mins_para = n_4mins_para/T_4mins_para
    power_4mins_para = (np.abs(X_4mins_para)**2)*sample_rate
    
    #find first number in x array higher than this limit and then last one lower
    #and then use that to find y-range to integrate over!
    x_4mins_where = np.where((freq_4mins_para > int_lower_lim) & (freq_4mins_para < int_upper_lim))
    
    x_4mins_toint = []
    y_4mins_para_toint = []
    y_4mins_perp_toint = []
    
    for i in x_4mins_where[0]:
        x_val = freq_4mins_para[i]
        y_para_val = power_4mins_para[i]
        y_perp_val = power_4mins_perp[i]
        x_4mins_toint.append(x_val)
        y_4mins_para_toint.append(y_para_val)
        y_4mins_perp_toint.append(y_perp_val)
    
    fourminute_para_int_power = np.trapz(y_4mins_para_toint, x_4mins_toint)
    fourminute_perp_int_power = np.trapz(y_4mins_perp_toint, x_4mins_toint)
    print('Power in 7-100mHz Pc3-4 Band, Four Minute Interval, Parallel:', fourminute_para_int_power, 'nT^2')
    print('Power in 7-100mHz Pc3-4 Band, Four Minute Interval, Perpendicular:', fourminute_perp_int_power, 'nT^2')
    
    #######TWO MINUTES
    
    #now, find average Cluster magnetic field direction during this time

    #Bx_gse, By_gse, Bz_gse, B_mag

    #normalised columns for each, then average

    ULF_df_2mins['Norm_Bx'] = ULF_df_2mins['Bx_gse'].div(ULF_df_2mins['B_mag'])
    ULF_df_2mins['Norm_By'] = ULF_df_2mins['By_gse'].div(ULF_df_2mins['B_mag'])
    ULF_df_2mins['Norm_Bz'] = ULF_df_2mins['Bz_gse'].div(ULF_df_2mins['B_mag'])

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
    Norm_2 = np.cross(B_mean, Norm_1)
    Norm_2 = Norm_2/(np.linalg.norm(Norm_2))
    #find B in direction of B perp 
    #first make array for each row in table, then get dot products for each 
    
    print('B mean 2mins:',   B_mean)
    print('Norm 1 2mins:',  Norm_1)
    print('Norm 2 2mins:',  Norm_2)
    
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
    
    para_mean_2 = ULF_df_2mins['B Para'].mean()
    ULF_df_2mins['B Para'] = ULF_df_2mins['B Para'] - para_mean_2

    #now calculate power spectrum for each perp component & sum.

    #FFT 2 mins norm 1
    x_2mins_perp_1 = ULF_df_2mins['B Perp 1'].to_numpy()

    X_2mins_perp_1 = rfft(x_2mins_perp_1,norm='ortho')
    N_2mins_perp_1 = len(X_2mins_perp_1)
    n_2mins_perp_1 = np.arange(N_2mins_perp_1)
    T_2mins_perp_1 = N_2mins_perp_1/sr
    freq_2mins_perp_1 = n_2mins_perp_1/T_2mins_perp_1
    power_2mins_perp_1 = (np.abs(X_2mins_perp_1)**2)*sample_rate

    #FFT 2 mins norm 2
    x_2mins_perp_2 = ULF_df_2mins['B Perp 2'].to_numpy()

    X_2mins_perp_2 = rfft(x_2mins_perp_2,norm='ortho')
    N_2mins_perp_2 = len(X_2mins_perp_2)
    n_2mins_perp_2 = np.arange(N_2mins_perp_2)
    T_2mins_perp_2 = N_2mins_perp_2/sr
    freq_2mins_perp_2 = n_2mins_perp_2/T_2mins_perp_2
    power_2mins_perp_2 = (np.abs(X_2mins_perp_2)**2)*sample_rate


    #add up (since power, do not need to sqrt). wait but freqs might not be same. 

    power_2mins_perp = power_2mins_perp_1 + power_2mins_perp_2

    #power_tot = power_2_1 + power_2_2
    
    #FFT parallel
    x_2mins_para = ULF_df_2mins['B Para'].to_numpy()

    X_2mins_para = rfft(x_2mins_para,norm='ortho')
    N_2mins_para = len(X_2mins_para)
    n_2mins_para = np.arange(N_2mins_para)
    T_2mins_para = N_2mins_para/sr
    freq_2mins_para = n_2mins_para/T_2mins_para
    power_2mins_para = (np.abs(X_2mins_para)**2)*sample_rate
    
    
    #lowfreq_lim = 1/120
    
    #add 5/3 line
    #starting at x = 10^-1 Hz and ending at 5.
    #starting y=30 
    A = 30/(0.1**(-5/3))
    x_array = np.geomspace(0.1, 5, num=50)
    y_array = A*(x_array**(-5/3))
    
    #find first number in x array higher than this limit and then last one lower
    #and then use that to find y-range to integrate over!
    x_2mins_where = np.where((freq_2mins_para > int_lower_lim) & (freq_2mins_para < int_upper_lim))
    
    x_2mins_toint = []
    y_2mins_para_toint = []
    y_2mins_perp_toint = []
    
    for i in x_2mins_where[0]:
        x_val = freq_2mins_para[i]
        y_para_val = power_2mins_para[i]
        y_perp_val = power_2mins_perp[i]
        x_2mins_toint.append(x_val)
        y_2mins_para_toint.append(y_para_val)
        y_2mins_perp_toint.append(y_perp_val)
    
    twominute_para_int_power = np.trapz(y_2mins_para_toint, x_2mins_toint)
    twominute_perp_int_power = np.trapz(y_2mins_perp_toint, x_2mins_toint)
    print('Power in 7-100mHz Pc3-4 Band, Two Minute Interval, Parallel:', twominute_para_int_power, 'nT^2')
    print('Power in 7-100mHz Pc3-4 Band, Two Minute Interval, Perpendicular:', twominute_perp_int_power, 'nT^2')
    
    #show spectra in two perp directions.
    fig, (ax3, ax4) = plt.subplots(2, figsize=(3,10))
    fig.suptitle('Interval Centred on ' + str_centre)
    
    ax3.set_title('4 Minute Interval')
    ax3.plot(freq_4mins_para, power_4mins_para, color="red", label='Parallel Power')
    ax3.plot(freq_4mins_perp_1, power_4mins_perp, color="black", label='Perpendicular Power')
    #ax3.plot(x_array, y_array, color="blue", label='-5/3 Slope')
    #plt.xlabel('Period (s)')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_ylabel('Power')
    ax3.set_xlabel('Frequency, Hz')
    ax3.set_xlim(0.001, 5)
    ax3.set_ylim(0.00001, 1000)
    #ax3.set_ylim(0.001, 100_000)
    ax3.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
    ax3.vlines(x=int_lower_lim, ymin = 0.0000001, ymax = 10_000, linestyles='dotted', color='k', label='Lower Bound')
    #ax3.vlines(x=dl_freq, ymin = 0.1, ymax = 100_000_000, linestyles='dashed', label='De Lauretis Freq')
    #ax3.vlines(x=tak_freq, ymin = 0.1, ymax = 100_000_000, linestyles='dotted', label='Takahashi Freq')
    #ax3.vlines(x=lowfreq_lim, ymin = 0.1, ymax = 100_000_000, linestyles='dotted', color='k', label='2 min window freq')
    ax3.vlines(x=int_upper_lim, ymin = 0.0000001, ymax = 10_000, linestyles='dotted', color='k', label='Upper Bound')
    #ax3.legend(loc='lower center')
    ax3.set_axisbelow(True)
    ax3.yaxis.grid(color='lightgray')
    ax3.xaxis.grid(color='lightgray')
    
    
    ax4.set_title('2 Minute Interval')
    ax4.plot(freq_2mins_para, power_2mins_para, color="red", label='Parallel Power')
    ax4.plot(freq_2mins_perp_1, power_2mins_perp, color="black", label='Perpendicular Power')
    #ax4.plot(x_array, y_array, color="blue", label='-5/3 Slope')
    #plt.xlabel('Period (s)')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_xlabel('Frequency, Hz')
    ax4.set_xlim(0.001, 5)
    ax4.set_ylim(0.00001, 1000)
    #ax4.set_ylim(0.001, 100_000)
    ax4.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
    ax4.vlines(x=int_lower_lim, ymin = 0.0000001, ymax = 10_000, linestyles='dotted', color='k', label='Lower Bound')
    #ax4.vlines(x=tak_freq, ymin = 0.1, ymax = 100_000_000, linestyles='dotted', label='Takahashi Freq')
    ax4.vlines(x=int_upper_lim, ymin = 0.0000001, ymax = 10_000, linestyles='dotted', color='k', label='Upper Bound')
    #ax4.legend(loc='lower center')
    ax4.set_axisbelow(True)
    ax4.yaxis.grid(color='lightgray')
    ax4.xaxis.grid(color='lightgray')
    
    plt.tight_layout()
    plt.show()
    
