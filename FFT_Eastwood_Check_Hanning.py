#Fourier Transform B Perp edition

import pandas as pd
import numpy as np
from numpy.fft import fft, ifft, rfft
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.ticker as ticker

pd.options.mode.chained_assignment = None

def FFT_Hann_20(cluster_ULF_csv, str_centre):
    
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
    time_start_10mins = time_cent - dt.timedelta(minutes=5)
    time_end_10mins = time_cent + dt.timedelta(minutes=5)    
    time_start_20mins = time_cent - dt.timedelta(minutes=10)
    time_end_20mins = time_cent + dt.timedelta(minutes=10)
    
    ULF_df_2mins = ULF_df.loc[((ULF_df.index >= time_start_2mins) & (ULF_df.index < time_end_2mins))]
    ULF_df_4mins = ULF_df.loc[((ULF_df.index >= time_start_4mins) & (ULF_df.index < time_end_4mins))]
    ULF_df_10mins = ULF_df.loc[((ULF_df.index >= time_start_10mins) & (ULF_df.index < time_end_10mins))]
    ULF_df_20mins = ULF_df.loc[((ULF_df.index >= time_start_20mins) & (ULF_df.index < time_end_20mins))]

    
    # sampling rate
    sr = 22
    sample_rate = 1/22
    ecf = np.sqrt(8/3)
    
    #######TWENTY MINUTES
    
    #now, find average Cluster magnetic field direction during this time

    #Bx_gse, By_gse, Bz_gse, B_mag

    mean_x = ULF_df_20mins['Bx_gse'].mean()
    mean_y = ULF_df_20mins['By_gse'].mean()
    mean_z = ULF_df_20mins['Bz_gse'].mean()

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
    
    list_of_B = []

    for i,j,k in zip(ULF_df_20mins['Bx_gse'], ULF_df_20mins['By_gse'], ULF_df_20mins['Bz_gse']):
        B_vec = np.array([i,j,k])
        list_of_B.append(B_vec)

    ULF_df_20mins['B Vector']= list_of_B

    list_norm_1_vals = []

    for i in ULF_df_20mins['B Vector']:
        projection_norm_1 = np.dot(i, Norm_1)
        list_norm_1_vals.append(projection_norm_1)

    ULF_df_20mins['B Perp 1'] = list_norm_1_vals

    list_norm_2_vals = []

    for i in ULF_df_20mins['B Vector']:
        projection_norm_2 = np.dot(i, Norm_2)
        list_norm_2_vals.append(projection_norm_2)
        
    list_para_vals = []

    for i in ULF_df_20mins['B Vector']:
        projection_para = np.dot(i, B_mean)
        list_para_vals.append(projection_para)

    ULF_df_20mins['B Perp 2'] = list_norm_2_vals
    ULF_df_20mins['B Para'] = list_para_vals
    
    para_mean = ULF_df_20mins['B Para'].mean()
    ULF_df_20mins['B Para'] = ULF_df_20mins['B Para'] - para_mean

    #now calculate power spectrum for each perp component & sum.
    #WITH HANN WINDOW
    sample_rate = 1/22
    ecf = np.sqrt(8/3)

    #FFT 2 mins norm 1
    x_20mins_perp_1 = ULF_df_20mins['B Perp 1'].to_numpy()
    twenty_min_hann = np.hanning(len(x_20mins_perp_1))
    x_20mins_perp_1 = x_20mins_perp_1*twenty_min_hann
    X_20mins_perp_1 = rfft(x_20mins_perp_1,norm='ortho')
    N_20mins_perp_1 = x_20mins_perp_1.size
    freq_20mins_perp_1 = np.fft.rfftfreq(N_20mins_perp_1, d=sample_rate)
    #ecf is energy correction for hann power loss
    power_20mins_perp_1 = 2*(np.abs(X_20mins_perp_1)**2)*ecf*sample_rate

    #FFT 2 mins norm 2
    x_20mins_perp_2 = ULF_df_20mins['B Perp 2'].to_numpy()
    x_20mins_perp_2 = x_20mins_perp_2*twenty_min_hann

    X_20mins_perp_2 = rfft(x_20mins_perp_2,norm='ortho')
    N_20mins_perp_2 = x_20mins_perp_2.size
    freq_20mins_perp_2 = np.fft.rfftfreq(N_20mins_perp_2, d=sample_rate)
    power_20mins_perp_2 = 2*(np.abs(X_20mins_perp_2)**2)*ecf*sample_rate


    #add up (since power, do not need to sqrt). wait but freqs might not be same. 

    power_20mins_perp = power_20mins_perp_1 + power_20mins_perp_2

    #power_tot = power_2_1 + power_2_2
    
    #FFT parallel
    x_20mins_para = ULF_df_20mins['B Para'].to_numpy()
    x_20mins_para = x_20mins_para*twenty_min_hann

    X_20mins_para = rfft(x_20mins_para,norm='ortho')
    N_20mins_para= x_20mins_para.size
    freq_20mins_para = np.fft.rfftfreq(N_20mins_para, d=sample_rate)
    power_20mins_para = 2*(np.abs(X_20mins_para)**2)*ecf*sample_rate
    
    
    int_lower_lim = 7*(10**(-3))
    int_upper_lim = 100*(10**(-3))
    
    #find first number in x array higher than this limit and then last one lower
    #and then use that to find y-range to integrate over!
    x_20mins_where = np.where((freq_20mins_para > int_lower_lim) & (freq_20mins_para < int_upper_lim))
    
    x_20mins_toint = []
    y_20mins_para_toint = []
    y_20mins_perp_toint = []
    
    for i in x_20mins_where[0]:
        x_val = freq_20mins_para[i]
        y_para_val = power_20mins_para[i]
        y_perp_val = power_20mins_perp[i]
        x_20mins_toint.append(x_val)
        y_20mins_para_toint.append(y_para_val)
        y_20mins_perp_toint.append(y_perp_val)
        
    twentyminute_para_int_power = np.trapz(y_20mins_para_toint, x_20mins_toint)
    twentyminute_perp_int_power = np.trapz(y_20mins_perp_toint, x_20mins_toint)

    #######TEN MINUTES
    
    #now, find average Cluster magnetic field direction during this time

    #Bx_gse, By_gse, Bz_gse, B_mag

    mean_x = ULF_df_10mins['Bx_gse'].mean()
    mean_y = ULF_df_10mins['By_gse'].mean()
    mean_z = ULF_df_10mins['Bz_gse'].mean()

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
    
    list_of_B = []

    for i,j,k in zip(ULF_df_10mins['Bx_gse'], ULF_df_10mins['By_gse'], ULF_df_10mins['Bz_gse']):
        B_vec = np.array([i,j,k])
        list_of_B.append(B_vec)

    ULF_df_10mins['B Vector']= list_of_B

    list_norm_1_vals = []

    for i in ULF_df_10mins['B Vector']:
        projection_norm_1 = np.dot(i, Norm_1)
        list_norm_1_vals.append(projection_norm_1)

    ULF_df_10mins['B Perp 1'] = list_norm_1_vals

    list_norm_2_vals = []

    for i in ULF_df_10mins['B Vector']:
        projection_norm_2 = np.dot(i, Norm_2)
        list_norm_2_vals.append(projection_norm_2)
        
    list_para_vals = []

    for i in ULF_df_10mins['B Vector']:
        projection_para = np.dot(i, B_mean)
        list_para_vals.append(projection_para)

    ULF_df_10mins['B Perp 2'] = list_norm_2_vals
    ULF_df_10mins['B Para'] = list_para_vals
    
    para_mean_10 = ULF_df_10mins['B Para'].mean()
    ULF_df_10mins['B Para'] = ULF_df_10mins['B Para'] - para_mean_10

    #now calculate power spectrum for each perp component & sum.

    #FFT 2 mins norm 1
    x_10mins_perp_1 = ULF_df_10mins['B Perp 1'].to_numpy()
    ten_min_hann = np.hanning(len(x_10mins_perp_1))
    x_10mins_perp_1 = x_10mins_perp_1*ten_min_hann

    X_10mins_perp_1 = rfft(x_10mins_perp_1,norm='ortho')
    N_10mins_perp_1= x_10mins_perp_1.size
    freq_10mins_perp_1 = np.fft.rfftfreq(N_10mins_perp_1, d=sample_rate)
    power_10mins_perp_1 = 2*(np.abs(X_10mins_perp_1)**2)*ecf*sample_rate

    #FFT 10 mins norm 2
    x_10mins_perp_2 = ULF_df_10mins['B Perp 2'].to_numpy()
    x_10mins_perp_2 = x_10mins_perp_2*ten_min_hann
    
    X_10mins_perp_2 = rfft(x_10mins_perp_2,norm='ortho')
    N_10mins_perp_2 = x_10mins_perp_2.size
    freq_10mins_perp_2 = np.fft.rfftfreq(N_10mins_perp_2, d=sample_rate)
    power_10mins_perp_2 = 2*(np.abs(X_10mins_perp_2)**2)*ecf*sample_rate

    #add up (since power, do not need to sqrt). wait but freqs might not be same. 

    power_10mins_perp = power_10mins_perp_1 + power_10mins_perp_2
    
    #FFT parallel
    x_10mins_para = ULF_df_10mins['B Para'].to_numpy()
    x_10mins_para = x_10mins_para*ten_min_hann
    X_10mins_para = rfft(x_10mins_para,norm='ortho')
    N_10mins_para= x_10mins_para.size
    freq_10mins_para = np.fft.rfftfreq(N_10mins_para, d=sample_rate)
    power_10mins_para = 2*(np.abs(X_10mins_para)**2)*ecf*sample_rate
    
    #find first number in x array higher than this limit and then last one lower
    #and then use that to find y-range to integrate over!
    x_10mins_where = np.where((freq_10mins_para > int_lower_lim) & (freq_10mins_para < int_upper_lim))
    
    x_10mins_toint = []
    y_10mins_para_toint = []
    y_10mins_perp_toint = []
    
    for i in x_10mins_where[0]:
        x_val = freq_10mins_para[i]
        y_para_val = power_10mins_para[i]
        y_perp_val = power_10mins_perp[i]
        x_10mins_toint.append(x_val)
        y_10mins_para_toint.append(y_para_val)
        y_10mins_perp_toint.append(y_perp_val)
    
    tenminute_para_int_power = np.trapz(y_10mins_para_toint, x_10mins_toint)
    tenminute_perp_int_power = np.trapz(y_10mins_perp_toint, x_10mins_toint)
    
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
    N_4mins_perp_1= x_4mins_perp_1.size
    freq_4mins_perp_1 = np.fft.rfftfreq(N_4mins_perp_1, d=sample_rate)
    power_4mins_perp_1 = 2*(np.abs(X_4mins_perp_1)**2)*ecf*sample_rate

    #FFT 4 mins norm 2
    x_4mins_perp_2 = ULF_df_4mins['B Perp 2'].to_numpy()
    x_4mins_perp_2 = x_4mins_perp_2*four_min_hann
    X_4mins_perp_2 = rfft(x_4mins_perp_2,norm='ortho')
    N_4mins_perp_2 = x_4mins_perp_2.size
    freq_4mins_perp_2 = np.fft.rfftfreq(N_4mins_perp_2, d=sample_rate)
    power_4mins_perp_2 = 2*(np.abs(X_4mins_perp_2)**2)*ecf*sample_rate


    #add up (since power, do not need to sqrt). wait but freqs might not be same. 

    power_4mins_perp = power_4mins_perp_1 + power_4mins_perp_2

    #power_tot = power_2_1 + power_2_2
    
    #FFT parallel
    x_4mins_para = ULF_df_4mins['B Para'].to_numpy()
    x_4mins_para = x_4mins_para*four_min_hann

    X_4mins_para = rfft(x_4mins_para,norm='ortho')
    N_4mins_para= x_4mins_para.size
    freq_4mins_para = np.fft.rfftfreq(N_4mins_para, d=sample_rate)
    power_4mins_para = 2*(np.abs(X_4mins_para)**2)*ecf*sample_rate

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
    
     #######TWO MINUTES
    
    #now, find average Cluster magnetic field direction during this time

    #Bx_gse, By_gse, Bz_gse
    
    mean_x = ULF_df_2mins['Bx_gse'].mean()
    mean_y = ULF_df_2mins['By_gse'].mean()
    mean_z = ULF_df_2mins['Bz_gse'].mean()

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
    two_min_hann = np.hanning(len(x_2mins_perp_1))
    x_2mins_perp_1 = x_2mins_perp_1*two_min_hann

    X_2mins_perp_1 = rfft(x_2mins_perp_1,norm='ortho')
    N_2mins_perp_1= x_2mins_perp_1.size
    freq_2mins_perp_1 = np.fft.rfftfreq(N_2mins_perp_1, d=sample_rate)
    power_2mins_perp_1 = 2*(np.abs(X_2mins_perp_1)**2)*ecf*sample_rate

    #FFT 2 mins norm 2
    x_2mins_perp_2 = ULF_df_2mins['B Perp 2'].to_numpy()
    x_2mins_perp_2 = x_2mins_perp_2*two_min_hann

    X_2mins_perp_2 = rfft(x_2mins_perp_2,norm='ortho')
    N_2mins_perp_2 = x_2mins_perp_2.size
    freq_2mins_perp_2 = np.fft.rfftfreq(N_2mins_perp_2, d=sample_rate)
    power_2mins_perp_2 = 2*(np.abs(X_2mins_perp_2)**2)*ecf*sample_rate


    #add up (since power, do not need to sqrt). wait but freqs might not be same. 

    power_2mins_perp = power_2mins_perp_1 + power_2mins_perp_2

    #power_tot = power_2_1 + power_2_2
    
    #FFT parallel
    x_2mins_para = ULF_df_2mins['B Para'].to_numpy()
    x_2mins_para = x_2mins_para*two_min_hann
    
    X_2mins_para = rfft(x_2mins_para,norm='ortho')
    N_2mins_para = x_2mins_para.size
    freq_2mins_para = np.fft.rfftfreq(N_2mins_para, d=sample_rate)
    power_2mins_para = 2*(np.abs(X_2mins_para)**2)*ecf*sample_rate
    
    #add 5/3 line
    #starting at x = 10^-1 Hz and ending at 5.
    #starting y=30 
    A = 30/(0.1**(-5/3))
    x_array = np.geomspace(0.1, 5, num=50)
    y_array = A*(x_array**(-5/3))
    
    #window = np.hanning(N_20mins_para)
    
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
    
    compressive_power_array = [twominute_para_int_power, fourminute_para_int_power, tenminute_para_int_power, twentyminute_para_int_power]
    transverse_power_array = [twominute_perp_int_power, fourminute_perp_int_power, tenminute_perp_int_power, twentyminute_perp_int_power]
    
    #show spectra in two perp directions.
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(5,16))
    fig.suptitle('Interval Centred on ' + str_centre)
    ax1.set_title('20 Minute Interval, Hann Window')
    ax1.plot(freq_20mins_para, power_20mins_para, color="red", label='Parallel Power')
    ax1.plot(freq_20mins_perp_1, power_20mins_perp, color="black", label='Perpendicular Power')
    #ax1.plot(x_array, y_array, color="blue", label='-5/3 Slope')
    #plt.xlabel('Period (s)')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel('Power')
    ax1.set_xlabel('Frequency, Hz')
    ax1.set_xlim(0.001, 5)
    ax1.set_ylim(0.0000001, 1000)
    ax1.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
    ax1.vlines(x=int_lower_lim, ymin = 0.0000001, ymax = 1000, linestyles='dotted', color='k', label='Lower Bound')
    #ax1.vlines(x=dl_freq, ymin = 0.1, ymax = 100_000_000, linestyles='dashed', label='De Lauretis Freq')
    #ax1.vlines(x=tak_freq, ymin = 0.1, ymax = 100_000_000, linestyles='dotted', label='Takahashi Freq')
    #ax1.vlines(x=lowfreq_lim, ymin = 0.1, ymax = 100_000_000, linestyles='dotted', color='k', label='2 min window freq')
    ax1.vlines(x=int_upper_lim, ymin = 0.0000001, ymax = 1000, linestyles='dotted', color='k', label='Upper Bound')
    ax1.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)
    ax1.set_axisbelow(True)
    ax1.yaxis.grid(color='lightgray')
    ax1.xaxis.grid(color='lightgray')
    
    ax2.set_title('10 Minute Interval, Hann Window')
    ax2.plot(freq_10mins_para, power_10mins_para, color="red", label='Parallel Power')
    ax2.plot(freq_10mins_perp_1, power_10mins_perp, color="black", label='Perpendicular Power')
    #ax2.plot(x_array, y_array, color="blue", label='-5/3 Slope')
    #plt.xlabel('Period (s)')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel('Power')
    ax2.set_xlabel('Frequency, Hz')
    ax2.set_xlim(0.001, 5)
    #ax2.set_ylim(0.001, 100_000)
    ax2.set_ylim(0.0000001, 1000)
    ax2.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
    ax2.vlines(x=int_lower_lim, ymin = 0.0000001, ymax = 1000, linestyles='dotted', color='k', label='Lower Bound')
    #ax2.vlines(x=dl_freq, ymin = 0.1, ymax = 100_000_000, linestyles='dashed', label='De Lauretis Freq')
    #ax2.vlines(x=tak_freq, ymin = 0.1, ymax = 100_000_000, linestyles='dotted', label='Takahashi Freq')
    #ax2.vlines(x=lowfreq_lim, ymin = 0.1, ymax = 100_000_000, linestyles='dotted', color='k', label='2 min window freq')
    ax2.vlines(x=int_upper_lim, ymin = 0.0000001, ymax = 1000, linestyles='dotted', color='k', label='Upper Bound')
    #ax2.legend(loc='lower center')
    ax2.set_axisbelow(True)
    ax2.yaxis.grid(color='lightgray')
    ax2.xaxis.grid(color='lightgray')
    
    
    ax3.set_title('4 Minute Interval, Hann Window')
    ax3.plot(freq_4mins_para, power_4mins_para, color="red", label='Parallel Power')
    ax3.plot(freq_4mins_perp_1, power_4mins_perp, color="black", label='Perpendicular Power')
    #ax3.plot(x_array, y_array, color="blue", label='-5/3 Slope')
    #plt.xlabel('Period (s)')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_ylabel('Power')
    ax3.set_xlabel('Frequency, Hz')
    ax3.set_xlim(0.001, 5)
    ax3.set_ylim(0.0000001, 1000)
    #ax3.set_ylim(0.001, 100_000)
    ax3.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
    ax3.vlines(x=int_lower_lim, ymin = 0.0000001, ymax = 1000, linestyles='dotted', color='k', label='Lower Bound')
    #ax3.vlines(x=dl_freq, ymin = 0.1, ymax = 100_000_000, linestyles='dashed', label='De Lauretis Freq')
    #ax3.vlines(x=tak_freq, ymin = 0.1, ymax = 100_000_000, linestyles='dotted', label='Takahashi Freq')
    #ax3.vlines(x=lowfreq_lim, ymin = 0.1, ymax = 100_000_000, linestyles='dotted', color='k', label='2 min window freq')
    ax3.vlines(x=int_upper_lim, ymin = 0.0000001, ymax = 1000, linestyles='dotted', color='k', label='Upper Bound')
    #ax3.legend(loc='lower center')
    ax3.set_axisbelow(True)
    ax3.yaxis.grid(color='lightgray')
    ax3.xaxis.grid(color='lightgray')
    
    
    ax4.set_title('2 Minute Interval, Hann Window')
    ax4.plot(freq_2mins_para, power_2mins_para, color="red", label='Parallel Power')
    ax4.plot(freq_2mins_perp_1, power_2mins_perp, color="black", label='Perpendicular Power')
    #ax4.plot(x_array, y_array, color="blue", label='-5/3 Slope')
    #plt.xlabel('Period (s)')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_xlabel('Frequency, Hz')
    ax4.set_xlim(0.001, 5)
    ax4.set_ylim(0.0000001, 1000)
    #ax4.set_ylim(0.001, 100_000)
    ax4.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
    ax4.vlines(x=int_lower_lim, ymin = 0.0000001, ymax = 1000, linestyles='dotted', color='k', label='Lower Bound')
    #ax4.vlines(x=tak_freq, ymin = 0.1, ymax = 100_000_000, linestyles='dotted', label='Takahashi Freq')
    ax4.vlines(x=int_upper_lim, ymin = 0.0000001, ymax = 1000, linestyles='dotted', color='k', label='Upper Bound')
    #ax4.legend(loc='lower center')
    ax4.set_axisbelow(True)
    ax4.yaxis.grid(color='lightgray')
    ax4.xaxis.grid(color='lightgray')
    
    plt.tight_layout()
    plt.show()
    
    return(compressive_power_array, transverse_power_array)
