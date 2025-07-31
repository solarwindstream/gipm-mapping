#averaging time window of ten minutes, centered around midpoint of Cluster interval (bearing in mind ref time is window start). Using 4 minute time window.
#this function takes 1 min OMNI readings and a list of time intervals and produces df of omni average B, V and n values
#including components
#for 10 minute window centered on midpoint of time interval 
import datetime as dt
import pandas as pd
import numpy as np
import math

#from https://stackoverflow.com/questions/11620914/how-do-i-remove-nan-values-from-a-numpy-array
#temporarily converting to dataframe in order to use pandas built-in dropna function
#values puts it back to numpy

def dropna(arr, *args, **kwarg):
    assert isinstance(arr, np.ndarray)
    dropped=pd.DataFrame(arr).dropna(*args, **kwarg).values
    if arr.ndim==1:
        dropped=dropped.flatten()
    return dropped

def omni_seg(om_df, only_full_windows):
    time_back = dt.timedelta(minutes=3)
    time_forward = dt.timedelta(minutes=7)

    #mean
    omni_B_ave = []
    omni_V_ave = []
    omni_N_ave = []
    omni_Bx_ave = []
    omni_By_ave = []
    omni_Bz_ave = []
    omni_Vx_ave = []
    omni_Vy_ave = []
    omni_Vz_ave = []
    omni_MA_ave = []
    omni_CA_ave = []
    omni_sc_dist_ave = []
        
    #median
    omni_B_med = []
    omni_V_med = []
    omni_N_med = []
    omni_Bx_med = []
    omni_By_med = []
    omni_Bz_med = []
    omni_Vx_med = []
    omni_Vy_med = []
    omni_Vz_med = []
    omni_MA_med = []
    omni_CA_med = []
    omni_sc_dist_med = []
    
    for i in only_full_windows:
    
        start_time = i - time_back
        end_time = i + time_forward

        mask = om_df.loc[(om_df.index >= start_time) & (om_df.index < end_time)]

        mean_B = mask.loc[:,'B_mag'].mean()
        mean_V = mask.loc[:,'V_gse'].mean()
        mean_N = mask.loc[:,'Np'].mean() 
        mean_Bx = mask.loc[:,'Bx_gse'].mean()
        mean_By = mask.loc[:,'By_gse'].mean()
        mean_Bz = mask.loc[:,'Bz_gse'].mean()
        mean_Vx = mask.loc[:,'Vx_gse'].mean()
        mean_Vy = mask.loc[:,'Vy_gse'].mean()
        mean_Vz = mask.loc[:,'Vz_gse'].mean()
        mean_MA = mask.loc[:,'M_A'].mean()
        mean_CA = mask.loc[:,'cone angle'].mean()
        #finding mean offset from Earth-Sun line
        y_dist = mask.loc[:,'Sc_y_gse'].to_numpy()
        z_dist = mask.loc[:,'Sc_z_gse'].to_numpy()
        dist_arr = (y_dist**2 + z_dist**2)**0.5
        dist_mean = np.mean(dist_arr)
        dist_median = np.median(dist_arr)
        
        
        median_B = mask.loc[:,'B_mag'].median()
        median_V = mask.loc[:,'V_gse'].median()
        median_N = mask.loc[:,'Np'].median() 
        median_Bx = mask.loc[:,'Bx_gse'].median()
        median_By = mask.loc[:,'By_gse'].median()
        median_Bz = mask.loc[:,'Bz_gse'].median()
        median_Vx = mask.loc[:,'Vx_gse'].median()
        median_Vy = mask.loc[:,'Vy_gse'].median()
        median_Vz = mask.loc[:,'Vz_gse'].median()
        median_MA = mask.loc[:,'M_A'].median()
        median_CA = mask.loc[:,'cone angle'].median()

        omni_B_ave.append(mean_B)
        omni_V_ave.append(mean_V)
        omni_N_ave.append(mean_N)
        omni_Bx_ave.append(mean_Bx)
        omni_By_ave.append(mean_By)
        omni_Bz_ave.append(mean_Bz)
        omni_Vx_ave.append(mean_Vx)
        omni_Vy_ave.append(mean_Vy)
        omni_Vz_ave.append(mean_Vz)
        omni_MA_ave.append(mean_MA)
        omni_CA_ave.append(mean_CA)
        omni_sc_dist_ave.append(dist_mean)
        
        omni_B_med.append(median_B)
        omni_V_med.append(median_V)
        omni_N_med.append(median_N)
        omni_Bx_med.append(median_Bx)
        omni_By_med.append(median_By)
        omni_Bz_med.append(median_Bz)
        omni_Vx_med.append(median_Vx)
        omni_Vy_med.append(median_Vy)
        omni_Vz_med.append(median_Vz)
        omni_MA_med.append(median_MA)
        omni_CA_med.append(median_CA)
        omni_sc_dist_med.append(dist_median)
        
        #now also determine maximum cone angle change within the window.
        Bx = mask.loc[:,'Bx_gse']
        By = mask.loc[:,'By_gse']
        Bz = mask.loc[:,'Bz_gse']
        
        #array of all vectors in this mask
        mag_field_vectors = np.transpose(np.array([Bx,By,Bz]))
        #normalise all vectors
        mag_field_vectors = mag_field_vectors/np.linalg.norm(mag_field_vectors, axis = 1, keepdims = True)
        #drop all NaN rows:
        
        mag_field_vectors = dropna(mag_field_vectors)
        
        #now for each row, take dot product with every other row and record maximum deviation between angles
        
        rot_angles = []

        for i in range(np.shape(mag_field_vectors)[0]):
            for j in range(np.shape(mag_field_vectors)[0]):
                dot_prod = np.dot(mag_field_vectors[i, :], mag_field_vectors[j, :])
                angle_in_deg = np.rad2deg(np.arccos(dot_prod))
                rot_angles.append(angle_in_deg)

        max_angle_deviation = max(rot_angles)

    om_averages = pd.DataFrame({'datetime': only_full_windows, 'Np (mean)': omni_N_ave, 'B_mag (mean)': omni_B_ave, 'V_gse (mean)': omni_V_ave, 'B_X_gse (mean)': omni_Bx_ave, 'B_Y_gse (mean)': omni_By_ave, 'B_Z_gse (mean)': omni_Bz_ave, 'V_X_gse (mean)': omni_Vx_ave, 'V_Y_gse (mean)': omni_Vy_ave, 'V_Z_gse (mean)': omni_Vz_ave, 'M_A (mean)':omni_MA_ave, 'cone angle (mean)':omni_CA_ave,'Distance from X line (mean)': omni_sc_dist_mean, 'Np (median)': omni_N_med, 'B_mag (median)': omni_B_med, 'V_gse (median)': omni_V_med, 'B_X_gse (median)': omni_Bx_med, 'B_Y_gse (median)': omni_By_med, 'B_Z_gse (median)': omni_Bz_med, 'V_X_gse (median)': omni_Vx_med, 'V_Y_gse (median)': omni_Vy_med, 'V_Z_gse (median)': omni_Vz_med, 'M_A (median)':omni_MA_med, 'cone angle (median)':omni_CA_med, 'Distance from X line (median)': omni_sc_dist_med, 'max IMF angle deviation': max_angle_deviation})

    return(om_averages)

