#now with 4 min files, find spectra for each

import pandas as pd
import numpy as np
from numpy.fft import rfft
import datetime as dt

from C1_Cluster_CDF_conv import C1_cdf_conv
from C2_Cluster_CDF_conv import C2_cdf_conv
from C3_Cluster_CDF_conv import C3_cdf_conv
from C4_Cluster_CDF_conv import C4_cdf_conv

pd.options.mode.chained_assignment = None

#input 22vps CDFs

cdf_list_path = '/data/scratch/apx059/23_Years_Data/CDFs/CDFs_'
year_n = '2023'
list_of_cdfs = cdf_list_path + year_n + '.csv'
list_files = pd.read_csv(list_of_cdfs, header=None)
list_files = list_files.rename(columns={0:'fnames'})
list_cdfs = list_files['fnames'].to_list()
year_path = '/data/scratch/apx059/23_Years_Data/CDFs/CDFs_' + year_n + '/'

df_list_c1 = []
df_list_c2 = []
df_list_c3 = []
df_list_c4 = []

for i in list_cdfs:
    if 'C1' in i:
        fpath = year_path + i
        df_c1 = C1_cdf_conv(fpath)
        a = df_c1.empty
        if not a:
            df_list_c1.append(df_c1)
    if 'C2' in i:
        fpath = year_path + i
        df_c2 = C2_cdf_conv(fpath)
        a = df_c2.empty
        if not a:
            df_list_c2.append(df_c2)
    if 'C3' in i:   
        fpath = year_path + i
        df_c3 = C3_cdf_conv(fpath)
        a = df_c3.empty
        if not a:
            df_list_c3.append(df_c3) 
    if 'C4' in i:   
        fpath = year_path + i
        df_c4 = C4_cdf_conv(fpath)
        a = df_c4.empty
        if not a:
            df_list_c4.append(df_c4)
            
#and now also the relevant GIPM CSVs

csv_list = r'/data/home/apx059/gipm-mapping/2023_4mins_list.csv'
gipm_4mins_path = r'/data/scratch/apx059/23_Years_Data/CSVs/GIPM_4mins/'
csv_files = pd.read_csv(csv_list, header=None)
csv_files = csv_files.rename(columns={0:'fnames'})
list_gipm_csvs = csv_files['fnames'].to_list()

#peel off spacecraft label and save for later 
imported_csv_list_1 = []
imported_csv_list_2 = []
imported_csv_list_3 = []
imported_csv_list_4 = []

for i in list_gipm_csvs:
    
    if 'C1' in i and not 'OMNI' in i:
        file_name = gipm_4mins_path + i
        df = pd.read_csv(file_name)
        df['datetime'] = pd.to_datetime(df['datetime'])
        df = df.set_index('datetime')
        imported_csv_list_1.append(df)
        
    if 'C2' in i and not 'OMNI' in i:
        file_name = gipm_4mins_path + i
        df = pd.read_csv(file_name)
        df['datetime'] = pd.to_datetime(df['datetime'])
        df = df.set_index('datetime')
        imported_csv_list_2.append(df)
        
    if 'C3' in i and not 'OMNI' in i:
        file_name = gipm_4mins_path + i
        df = pd.read_csv(file_name)
        df['datetime'] = pd.to_datetime(df['datetime'])
        df = df.set_index('datetime')
        imported_csv_list_3.append(df)
        
    if 'C4' in i and not 'OMNI' in i:
        file_name = gipm_4mins_path + i
        df = pd.read_csv(file_name)
        df['datetime'] = pd.to_datetime(df['datetime'])
        df = df.set_index('datetime')
        imported_csv_list_4.append(df)

#now for each spacecraft, concatenate all the dataframes (gipm and raw separately), and then run through the datetimes in the GIPM df
#running FFT for each, using raw B data from cdfs

#first define function for one reading

def FFT_Hann(cluster_raw_df, window_start):
    
    #mask to just relevant times

    time_start_4mins = window_start
    time_end_4mins = window_start + dt.timedelta(minutes=4)

    ULF_df_4mins = cluster_raw_df.loc[((cluster_raw_df.index >= time_start_4mins) & (cluster_raw_df.index < time_end_4mins))]

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

    for i in ULF_df_4mins['B Vector']
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

##############
#now a function that will run FFT_Hann 
#the index should be the datetime

def FFT_run(gipm_df, raw_df, sc_string):
    
    para_i_p_list = []
    perp_i_p_list = []
    
    for i in gipm_df.index:
        para_i_p, perp_i_p, freq, p_para, p_perp_1, p_perp_2 = FFT_Hann(raw_df, i)
        para_i_p_list.append(para_i_p)
        perp_i_p_list.append(perp_i_p)
        spectral_df = pd.DataFrame({'Freq':freq, 'Parallel Power':p_para, 'Perp 1 Power': p_perp_1, 'Perp 2 Power':p_perp_2})
        df_name = str(i) + sc_string 
        df_path = '/data/scratch/apx059/23_Years_Data/CSVs/GIPM_4mins/Fourier_Products/'
        fpath = df_path + df_name
        spectral_df.to_csv(fpath)
     
    #save the perp and para integrated values. add on for each window in original csv?
    para_array = np.array(para_i_p_list)
    perp_array = np.array(perp_i_p_list)
    datetime_list = gipm_df.index
    power_df = pd.DataFrame({'datetime': datetime_list, 'Integrated Parallel Power':para_array, 'Integrated Perpendicular Power':perp_array})
    power_df['datetime'] = pd.to_datetime(power_df['datetime'])
    power_df = power_df.set_index('datetime')
    
    return power_df
    
    
    
#now run for each SC

gipm_df_1 = pd.concat(imported_csv_list_1)
raw_df_1 = pd.concat(df_list_c1)
gipm_df_2 = pd.concat(imported_csv_list_2)
raw_df_2 = pd.concat(df_list_c2)
gipm_df_3 = pd.concat(imported_csv_list_3)
raw_df_3 = pd.concat(df_list_c3)
gipm_df_4 = pd.concat(imported_csv_list_4)
raw_df_4 = pd.concat(df_list_c4)

power_1 = FFT_run(gipm_df_1,raw_df_1,'C1')
power_2 = FFT_run(gipm_df_2,raw_df_2,'C2')
power_3 = FFT_run(gipm_df_3,raw_df_3,'C3')
power_4 = FFT_run(gipm_df_4,raw_df_4,'C4')

gipm_df_up_1 = pd.merge(gipm_df_1,power_1, left_index=True, right_index=True)
gipm_df_up_2 = pd.merge(gipm_df_2,power_2, left_index=True, right_index=True)
gipm_df_up_3 = pd.merge(gipm_df_3,power_3, left_index=True, right_index=True)
gipm_df_up_4 = pd.merge(gipm_df_4,power_4, left_index=True, right_index=True)

fpath_1 = '/data/scratch/apx059/23_Years_Data/CSVs/GIPM_4mins/Fourier_Products/2023_C1.csv'
fpath_2 = '/data/scratch/apx059/23_Years_Data/CSVs/GIPM_4mins/Fourier_Products/2023_C2.csv'
fpath_3 = '/data/scratch/apx059/23_Years_Data/CSVs/GIPM_4mins/Fourier_Products/2023_C3.csv'
fpath_4 = '/data/scratch/apx059/23_Years_Data/CSVs/GIPM_4mins/Fourier_Products/2023_C4.csv'

gipm_df_up_1.to_csv(fpath_1)
gipm_df_up_2.to_csv(fpath_2)
gipm_df_up_3.to_csv(fpath_3)
gipm_df_up_4.to_csv(fpath_4)

