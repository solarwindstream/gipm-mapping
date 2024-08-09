#for every segment, transform location using given coeffs, for later plotting. Just C1 for now
import numpy as np
from new_xyz import new_xyz
import datetime as dt
import pandas as pd

def gipm_loc_transform(only_full_windows, df_a, GIPM_matrices, FAC_coeffs):
    
    time_window = dt.timedelta(seconds=120)
    Cluster_GIPM_orbits = []
    
    #mask to just the two minute window in cluster data and find mean location
    
    for i,k in zip(only_full_windows, FAC_coeffs):
        start_time = i
        end_time = i + time_window
        mask = df_a.loc[(df_a.index >= start_time) & (df_a.index < end_time)]
        m = FAC_coeffs.index(k)
        
        x_gse_mean = mask['X_gse'].mean()
        y_gse_mean = mask['Y_gse'].mean()
        z_gse_mean = mask['Z_gse'].mean()
        
        #form location arrays for every observation and then transform using coeffs/FAC
        Cluster_GSE = np.array([x_gse_mean, y_gse_mean, z_gse_mean])
        #rotate
        trans_mat = GIPM_matrices[m]
        Cluster_GIPM = new_xyz(Cluster_GSE, trans_mat)
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
                
        Cluster_GIPM = np.array([x_a, y_a, z_a])
                
        #append result
        Cluster_GIPM_orbits.append(Cluster_GIPM)
        
    Cluster_dt_loc = pd.DataFrame({'datetime':only_full_windows, 'GIPM Loc': Cluster_GIPM_orbits})
        

    return(Cluster_dt_loc)  
    
