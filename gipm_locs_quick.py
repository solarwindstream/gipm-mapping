#now transform all locations to GIPM
#input cluster GSE-location dataframes matched with GIPM matrices & FAC Coefficients
#output list of cluster GIPM-location
import pandas as pd
import numpy as np
from new_xyz import new_xyz

def gipm_locs_quick(df_cluster, GIPM_matrices, FAC_coeffs):
    
    Cluster_GIPM_orbits = []
    df_cluster_ind = df_cluster.index
    #mask to just the two minute window in cluster data and find mean location
    
    for i,k in zip(df_cluster_ind, FAC_coeffs):
        m = FAC_coeffs.index(k)
        
        x_gse = df_cluster.loc[i,'X_gse']
        y_gse = df_cluster.loc[i,'Y_gse']
        z_gse = df_cluster.loc[i,'Z_gse']
        
        #form location arrays for every observation and then transform using coeffs/FAC
        Cluster_GSE = np.array([x_gse, y_gse, z_gse])
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
        
    Cluster_dt_loc = pd.DataFrame({'datetime':i, 'GIPM Loc': Cluster_GIPM_orbits})
        

    return(Cluster_dt_loc)  
    

    
