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

    #Bx_gse, By_gse, Bz_gse

    mean_x = ULF_df_4mins['Bx_gse'].mean()
    mean_y = ULF_df_4mins['By_gse'].mean()
    mean_z = ULF_df_4mins['Bz_gse'].mean()

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
    N_4mins_perp_1= x_4mins_perp_1.size
    freq_4mins_perp_1 = np.fft.rfftfreq(N_4mins_perp_1, d=sample_rate)
    power_4mins_perp_1 = 2*(np.abs(X_4mins_perp_1)**2)*ecf*sample_rate

    #FFT 2 mins norm 2
    x_4mins_perp_2 = ULF_df_4mins['B Perp 2'].to_numpy()
    x_4mins_perp_2 = x_4mins_perp_2*four_min_hann

    X_4mins_perp_2 = rfft(x_4mins_perp_2,norm='ortho')
    N_4mins_perp_2 = x_4mins_perp_2.size
    freq_4mins_perp_2 =np.fft.rfftfreq(N_4mins_perp_2,d=sample_rate)
    power_4mins_perp_2 = 2*(np.abs(X_4mins_perp_2)**2)*ecf*sample_rate

    #add up (since power, do not need to sqrt). wait but freqs might not be same. 

    power_4mins_perp = power_4mins_perp_1 + power_4mins_perp_2

    #power_tot = power_2_1 + power_2_2
    
    #FFT parallel
    x_4mins_para = ULF_df_4mins['B Para'].to_numpy()
    x_4mins_para = x_4mins_para*four_min_hann

    X_4mins_para = rfft(x_4mins_para,norm='ortho')
    N_4mins_para = x_4mins_para.size
    freq_4mins_para =np.fft.rfftfreq(N_4mins_para,d=sample_rate)
    power_4mins_para = 2*(np.abs(X_4mins_para)**2)*ecf*sample_rate
    
    ## Compute the power into the required frequency range
    F_LOW,F_HIGH=7e-3,1e-1
    mask=(freq_4mins_para>F_LOW) & (freq_4mins_para<F_HIGH)
    delta_f=freq_4mins_para[1]-freq_4mins_para[0]
    P_para=np.sum(power_4mins_para[mask])*delta_f
    print('int para power:', P_para,'nT^2')
    
    F_LOW,F_HIGH=7e-3,1e-1
    mask=(freq_4mins_perp_1>F_LOW) & (freq_4mins_perp_1<F_HIGH)
    delta_f=freq_4mins_perp_1[1]-freq_4mins_perp_1[0]
    P_perp=np.sum(power_4mins_perp[mask])*delta_f
    print('int perp power:', P_perp,'nT^2')
    
    #show windowed vs unwindowed components
    
    fig, (ax0, ax1, ax2) = plt.subplots(3, figsize=(8, 12), tight_layout=True)
    
    ax0.set_title('Original Time Series')
    ax0.plot(ULF_df_4mins.index, ULF_df_4mins['Bx_gse'], color = 'g', label='$B_{x}$')
    ax0.plot(ULF_df_4mins.index, ULF_df_4mins['By_gse'], color = 'b', label='$B_{y}$')
    ax0.plot(ULF_df_4mins.index, ULF_df_4mins['Bz_gse'], color = 'r', label='$B_{z}$')
    ax0.plot(ULF_df_4mins.index, ULF_df_4mins['B_mag'], color = 'k', label='$B_{mag}$')
    ax0.set_xlabel('Time')
    ax0.set_ylabel('Magnetic Field, nT')
    ax0.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    ax0.set_axisbelow(True)
    ax0.xaxis.set_major_formatter(formatter)
    ax0.xaxis.set_minor_locator(AutoMinorLocator())
    ax0.yaxis.grid(color='lightgray')
    ax0.xaxis.grid(color='lightgray')
    
    #fig.suptitle('Interval Centred on ' + str_centre)
    ax1.set_title('Mean Removed')
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
    
    return(freq_4mins_para, power_4mins_para, power_4mins_perp)