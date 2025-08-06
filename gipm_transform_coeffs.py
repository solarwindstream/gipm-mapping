#gipm transform 
import numpy as np

def gipm_transform_coeffs_mean(om_averages):
    
    #find average MA for each whole interval first
    from math import pi
    from GIPM_Trans import GIPM_trans

    GIPM_X_Vecs = []
    GIPM_Y_Vecs = []
    GIPM_Z_Vecs = []
    FAC_coeffs = []

    #taking average B and V vectors from input, apply GIPM_trans function to find GIPM basis vectors
    #and scaling coefficient (FAC) during this time window.
    
    for i in om_averages.index:
        
        b_x = om_averages.loc[i,'B_X_gse (mean)']
        b_y = om_averages.loc[i,'B_Y_gse (mean)']
        b_z = om_averages.loc[i,'B_Z_gse (mean)']
        v_x = om_averages.loc[i,'V_X_gse (mean)']
        v_y = om_averages.loc[i,'V_Y_gse (mean)']
        v_z = om_averages.loc[i,'V_Z_gse (mean)']        
        b_gse = np.array([b_x,b_y,b_z])
        v_gse = np.array([v_x,v_y,v_z])
        X_GIPM, Y_GIPM, Z_GIPM = GIPM_trans(b_gse, v_gse)
        
        GIPM_X_Vecs.append(X_GIPM)
        GIPM_Y_Vecs.append(Y_GIPM)
        GIPM_Z_Vecs.append(Z_GIPM)
        
    #now for each element find FAC
    
    for i in om_averages.index:
        n_p = om_averages.loc[i,'Np (mean)']
        v_x = om_averages.loc[i,'V_X_gse (mean)']
        v_y = om_averages.loc[i,'V_Y_gse (mean)']
        v_z = om_averages.loc[i,'V_Z_gse (mean)']        
        v_gse = np.array([v_x,v_y,v_z])
        Vabs = np.linalg.norm(v_gse)
        FAC=((n_p/7.0)*(Vabs/457.5)**2)**0.166667
        
        FAC_coeffs.append(FAC)
        
    return(GIPM_X_Vecs, GIPM_Y_Vecs, GIPM_Z_Vecs, FAC_coeffs)

def gipm_transform_coeffs_median(om_averages):
    
    #find average MA for each whole interval first
    from math import pi
    from GIPM_Trans import GIPM_trans
    
    GIPM_X_Vecs = []
    GIPM_Y_Vecs = []
    GIPM_Z_Vecs = []
    FAC_coeffs =[]

    #now for each element in om_averages array, generate GIPM transform to transform Cluster points
    for i in om_averages.index:
        
        b_x = om_averages.loc[i,'B_X_gse (median)']
        b_y = om_averages.loc[i,'B_Y_gse (median)']
        b_z = om_averages.loc[i,'B_Z_gse (median)']
        v_x = om_averages.loc[i,'V_X_gse (median)']
        v_y = om_averages.loc[i,'V_Y_gse (median)']
        v_z = om_averages.loc[i,'V_Z_gse (median)']        
        b_gse = np.array([b_x,b_y,b_z])
        v_gse = np.array([v_x,v_y,v_z])
        X_GIPM, Y_GIPM, Z_GIPM = GIPM_trans(b_gse, v_gse)
        
        GIPM_X_Vecs.append(X_GIPM)
        GIPM_Y_Vecs.append(Y_GIPM)
        GIPM_Z_Vecs.append(Z_GIPM)
        
    #now for each element find FAC
    
    for i in om_averages.index:
        n_p = om_averages.loc[i,'Np (median)']
        v_x = om_averages.loc[i,'V_X_gse (median)']
        v_y = om_averages.loc[i,'V_Y_gse (median)']
        v_z = om_averages.loc[i,'V_Z_gse (median)']        
        v_gse = np.array([v_x,v_y,v_z])
        Vabs = np.linalg.norm(v_gse)
        FAC=((n_p/7.0)*(Vabs/457.5)**2)**0.166667
        
        FAC_coeffs.append(FAC)
        
    return(GIPM_X_Vecs, GIPM_Y_Vecs, GIPM_Z_Vecs, FAC_coeffs)
