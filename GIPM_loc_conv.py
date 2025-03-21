#for every segment, transform location using given coeffs, for later plotting.
import numpy as np
from new_xyz import new_xyz
import datetime as dt
import pandas as pd

def gipm_loc_transform(only_full_windows, df_a, GIPM_X_Vecs, GIPM_Y_Vecs, GIPM_Z_Vecs, FAC_coeffs):
    
    time_window = dt.timedelta(seconds=120)
    Cluster_GIPM_X = []
    Cluster_GIPM_Y = []
    Cluster_GIPM_Z = []
    
    #mask to just the two minute window in cluster data and find mean location
    
    for i,k,l,m,n in zip(only_full_windows, FAC_coeffs, GIPM_X_Vecs, GIPM_Y_Vecs, GIPM_Z_Vecs):
        start_time = i
        end_time = i + time_window
        mask = df_a.loc[(df_a.index >= start_time) & (df_a.index < end_time)]
        
        x_gse_mean = mask['X_gse'].mean()
        y_gse_mean = mask['Y_gse'].mean()
        z_gse_mean = mask['Z_gse'].mean()
        
        #form location arrays for every observation and then transform using coeffs/FAC
        Cluster_GSE = np.array([x_gse_mean, y_gse_mean, z_gse_mean])
        
        #print('X, Y, Z Vecs', l,m,n,'Cluster Loc', Cluster_GSE)
        
        #coefficients in new frame 
        
        X_coeff = np.dot(Cluster_GSE,l)
        Y_coeff = np.dot(Cluster_GSE,m)
        Z_coeff = np.dot(Cluster_GSE,n)
        
        #print('X_coeff shape', X_coeff.shape)
        #print('Y_coeff shape', Y_coeff.shape)
        #print('Z_coeff shape', Z_coeff.shape)
       
        Cluster_GIPM = np.array([X_coeff, Y_coeff, Z_coeff])
        
        #FAC scaling
        #r is scaled by this but the direction should remain the same
        x_a = Cluster_GIPM[0]
        y_a = Cluster_GIPM[1]
        z_a = Cluster_GIPM[2]
        r_a = (x_a**2 + y_a**2 + z_a**2)**0.5
        #find original angles
        theta_a = np.arccos(z_a/r_a)
        phi_a = np.arctan(y_a/x_a)

        r_a = r_a*k

        #new x_a, y_a, z_a
        x_a = r_a*(np.sin(theta_a))*(np.cos(phi_a))
        y_a = r_a*(np.sin(theta_a))*(np.sin(phi_a))
        z_a = r_a*(np.cos(theta_a))
                
        #Cluster_GIPM = np.array([x_a, y_a, z_a])
                
        #append result
        Cluster_GIPM_X.append(x_a)
        Cluster_GIPM_Y.append(y_a)
        Cluster_GIPM_Z.append(z_a)
        
    Cluster_dt_loc = pd.DataFrame({'datetime':only_full_windows, 'GIPM X': Cluster_GIPM_X, 'GIPM Y': Cluster_GIPM_Y, 'GIPM Z': Cluster_GIPM_Z})
        

    return(Cluster_dt_loc)  
    
