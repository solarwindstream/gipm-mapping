#gipm transform 
import numpy as np

#b_gse = np.array([B_X_ave,B_Y_ave,B_Z_ave])
#v_gse = np.array([V_x_mean,V_y_mean,V_z_mean])
#removed only_full_windows, df_a from inputs (not needed)

def gipm_transform(om_averages):
    #find average MA for each whole interval first
    from math import pi
    from GIPM_Trans import GIPM_trans
    #from merka05_surface_eq_array_GIPM_2 import merka05_surface_eq_array_GIPM

    #mu0 = 4*pi*10**(-7) #permittivity vs/(Am)
    #m_p = 1.672631*10**(-27) #proton mass kg

    #Babs = om_averages.loc[:,'B_mag'].mean()
    #Vabs = om_averages.loc[:,'V_gse'].mean()
    #Nave = om_averages.loc[:,'Np'].mean()
    #rho = Nave*10**6 #in #/m^3

    #VA = Babs*10**(-9)/np.sqrt(mu0*rho*m_p)/1000 #km/s
    #XMA = Vabs/VA
    
    #get bow shock shape for this average XMA
    #fitting_coeffs = merka05_surface_eq_array_GIPM(XMA)
    
    GIPM_matrices = []
    FAC_coeffs =[]

    #now for each element in om_averages array, generate GIPM transform to transform Cluster points
    for i in om_averages.index:
        
        b_x = om_averages.loc[i,'B_X_gse']
        b_y = om_averages.loc[i,'B_Y_gse']
        b_z = om_averages.loc[i,'B_Z_gse']
        v_x = om_averages.loc[i,'V_X_gse']
        v_y = om_averages.loc[i,'V_Y_gse']
        v_z = om_averages.loc[i,'V_Z_gse']        
        b_gse = np.array([b_x,b_y,b_z])
        v_gse = np.array([v_x,v_y,v_z])
        GIPM_t = GIPM_trans(b_gse, v_gse)
        
        GIPM_matrices.append(GIPM_t)
    
    #now for each element find FAC
    
    for i in om_averages.index:
        n_p = om_averages.loc[i,'Np']
        v_x = om_averages.loc[i,'V_X_gse']
        v_y = om_averages.loc[i,'V_Y_gse']
        v_z = om_averages.loc[i,'V_Z_gse']        
        v_gse = np.array([v_x,v_y,v_z])
        Vabs = np.linalg.norm(v_gse)
        FAC=((n_p/7.0)*(Vabs/457.5)**2)**0.166667
        
        FAC_coeffs.append(FAC)
        
    return(GIPM_matrices, FAC_coeffs)
