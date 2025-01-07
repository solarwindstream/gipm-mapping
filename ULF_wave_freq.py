##code to determine expected ulf wave frequency for proton
import math
import numpy as np
from merka05_closest_point_GIPM import merka05_closest_point
from GIPM_Trans import GIPM_trans
import pandas as pd

def ULF_freq(XMA, cluster_gipm, om_df, int_start):

    m_p =  1.67262192 * 10**(-27)
    q = 1.60217663 * 10**(-19)

    #nb this requires om_df only containing one copy of each instance. don't stick them together.
    slice_om = om_df.loc[om_df.index == int_start]
    B_nT = slice_om.loc[int_start,'B_mag']
    #B = B_nT*(10**(-9))

    pi = np.pi
    pref = (q*(10**(-6)))/(2*m_p*2*pi)
    
    B_x = om_df.loc[int_start,'Norm Bx']
    B_y = om_df.loc[int_start,'Norm By']
    B_z = om_df.loc[int_start,'Norm Bz']
    
    B_vec = np.array([B_x, B_y, B_z])
    B_vec = B_vec/(np.linalg.norm(B_vec))
    
    #directions
    theta_B_x = om_df.loc[int_start,'cone angle']
    print(theta_B_x)
    cos_theta_B_x = np.cos(np.deg2rad(theta_B_x))

    r_bs_gipm, normal_gipm = merka05_closest_point(XMA, cluster_gipm)
    #this is in GIPM but should we be re-adjusting for GSE?
    cos_theta_n_X = np.dot(normal_gipm, [1,0,0])
    print('cos_theta_nX:', cos_theta_n_X)
    
    #doing this in GIPM frame --> rotate B into GIPM frame
    
    b_gse_x = slice_om.loc[int_start,'B_X_gse']
    b_gse_y = slice_om.loc[int_start,'B_Y_gse']
    b_gse_z = slice_om.loc[int_start,'B_Z_gse']
    v_gse_x = slice_om.loc[int_start,'V_X_gse']
    v_gse_y = slice_om.loc[int_start,'V_Y_gse']
    v_gse_z = slice_om.loc[int_start,'V_Z_gse']
    
    b_gse = np.array([b_gse_x, b_gse_y, b_gse_z])
    v_gse = np.array([v_gse_x, v_gse_y, v_gse_z])
    
    X_GIPM, Y_GIPM, Z_GIPM = GIPM_trans(b_gse, v_gse)
    
    b_gipm_x = np.dot(b_gse, X_GIPM)
    b_gipm_y = np.dot(b_gse, Y_GIPM)
    b_gipm_z = np.dot(b_gse, Z_GIPM)
    
    b_gipm = np.array([b_gipm_x, b_gipm_y, b_gipm_z])
    b_gipm = b_gipm/np.linalg.norm(b_gipm)
    
    cos_theta_B_n = np.dot(b_gipm, normal_gipm)
    print('cos_theta_B_n:', cos_theta_B_n)

    #omega is pref x cos(thetaB_X) x cos(thetaB_n)/cos(thetan_x)
    
    #returns results in mHz
    dl_freq = (pref*B_nT*cos_theta_B_x*cos_theta_B_n)/cos_theta_n_X
    dl_freq = np.abs(dl_freq)
    
    tak_freq = (pref*B_nT*(cos_theta_B_x**2))
    tak_freq = np.abs(tak_freq)
    
    return dl_freq, tak_freq





