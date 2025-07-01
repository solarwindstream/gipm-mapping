import pandas as pd
import numpy as np
from numpy.fft import rfft
import datetime as dt

##Produces a four-minute power spectra in perp & parallel directions, using the Hann window. Also returns integrated power in the
##ULF wave band (7-100mHz)
#requires that the input is already masked to 4 mins!

def FFT_Hann(ULF_df_4mins):
    
    #mask to just relevant times

    #time_start_4mins = window_start
    #time_end_4mins = window_start + dt.timedelta(minutes=4)

    #ULF_df_4mins = cluster_raw_df.loc[((cluster_raw_df.index >= time_start_4mins) & (cluster_raw_df.index < time_end_4mins))]

    # sampling rate
    sr = 22
    sample_rate = 1/22
    ecf = np.sqrt(8/3)
    
    #pc3-4 pulsation window
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

    B_mean = np.array([mean_x, mean_y, mean_z])
    X_dir = np.array([1,0,0])

    Norm_1 = np.cross(B_mean, X_dir)
    Norm_1 = Norm_1/(np.linalg.norm(Norm_1))
    Norm_2 = np.cross(B_mean, Norm_1)
    Norm_2 = Norm_2/(np.linalg.norm(Norm_2))
    #find B in direction of B perp 
    #first make array for each row in table, then get dot products for each 

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

    #FFT 4 mins norm 1
    x_4mins_perp_1 = ULF_df_4mins['B Perp 1'].to_numpy()
    four_min_hann = np.hanning(len(x_4mins_perp_1))
    x_4mins_perp_1 = x_4mins_perp_1*four_min_hann

    X_4mins_perp_1 = rfft(x_4mins_perp_1,norm='ortho')
    N_4mins_perp_1 = len(X_4mins_perp_1)
    n_4mins_perp_1 = np.arange(N_4mins_perp_1)
    T_4mins_perp_1 = N_4mins_perp_1/sr
    freq_4mins_perp_1 = n_4mins_perp_1/T_4mins_perp_1
    power_4mins_perp_1 = (np.abs(X_4mins_perp_1)**2)*ecf*sample_rate

    #FFT 4 mins norm 2
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

    #FFT parallel
    x_4mins_para = ULF_df_4mins['B Para'].to_numpy()
    x_4mins_para = x_4mins_para*four_min_hann

    X_4mins_para = rfft(x_4mins_para,norm='ortho')
    N_4mins_para = len(X_4mins_para)
    n_4mins_para = np.arange(N_4mins_para)
    T_4mins_para = N_4mins_para/sr
    freq_4mins_para = n_4mins_para/T_4mins_para
    power_4mins_para = (np.abs(X_4mins_para)**2)*ecf*sample_rate

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

    fourminute_para_int_power = np.trapezoid(y_4mins_para_toint, x_4mins_toint)
    fourminute_perp_int_power = np.trapezoid(y_4mins_perp_toint, x_4mins_toint)
    
    return(fourminute_para_int_power, fourminute_perp_int_power, freq_4mins_para, power_4mins_para, power_4mins_perp_1, power_4mins_perp_2)