import cdflib
import pandas as pd
import numpy as np
import math
##read in OMNI data
#calculate cone angle
#produce OMNI dataframe with datetime index

def OMNI_cdf_conv(filename):
    cdf_file = cdflib.CDF(filename)

    datetimes_omni = cdflib.cdfepoch.encode(cdf_file['Epoch'])

    omni_dt_series = pd.Series(datetimes_omni) 
    datetimes_omni_a = pd.to_datetime(omni_dt_series)
    datetimes_omni_a = datetimes_omni_a.to_frame(name="datetime")

    om_df = pd.DataFrame({'Np': cdf_file['proton_density'], 'Bx_gse': cdf_file['BX_GSE'], 'By_gse': cdf_file['BY_GSE'], "Bz_gse": cdf_file['BZ_GSE'], 'B_mag': cdf_file['F'], "Vx_gse": cdf_file['Vx'], "Vy_gse": cdf_file['Vy'],"Vz_gse": cdf_file['Vz'], 'V_gse': cdf_file['flow_speed'], 'M_A': cdf_file['Mach_num'], 'Sc_x_gse': cdf_file['x'], 'Sc_y_gse': cdf_file['y'], 'Sc_z_gse': cdf_file['z'], 'BS_x_gse': cdf_file['BSN_x'], 'BS_y_gse': cdf_file['BSN_y'], 'BS_z_gse': cdf_file['BSN_z']})
    
    om_df = om_df.join(datetimes_omni_a)
    om_df = om_df.set_index('datetime')

    #replace fill vals
    om_df[om_df['Bx_gse'] == 9999.99] = np.nan
    om_df[om_df['By_gse'] == 9999.99] = np.nan
    om_df[om_df['Bz_gse'] == 9999.99] = np.nan
    om_df[om_df['B_mag'] == 9999.99] = np.nan
    om_df[om_df['Vx_gse'] == 99999.9] = np.nan
    om_df[om_df['Vy_gse'] == 99999.99] = np.nan
    om_df[om_df['Vz_gse'] == 99999.99] = np.nan
    om_df[om_df['V_gse'] == 99999.99] = np.nan
    om_df[om_df['M_A'] == 999.9] = np.nan
    om_df[om_df['Sc_x_gse'] == 9999.99] = np.nan
    om_df[om_df['Sc_y_gse'] == 9999.99] = np.nan
    om_df[om_df['Sc_z_gse'] == 9999.99] = np.nan
    om_df[om_df['BS_x_gse'] == 9999.99] = np.nan
    om_df[om_df['BS_y_gse'] == 9999.99] = np.nan
    om_df[om_df['BS_z_gse'] == 9999.99] = np.nan
    
    
    #now cone angle!
    omni_cone_angle = []
    error_list_NaN = []
    x_vec = np.array([1, 0, 0])
    
    for i in om_df.index:

        bx_temp = om_df.loc[i, 'Bx_gse']
        by_temp = om_df.loc[i, 'By_gse']
        bz_temp = om_df.loc[i, 'Bz_gse']
        bmag_temp = om_df.loc[i, 'B_mag']

        bx_normed = bx_temp/bmag_temp
        by_normed = by_temp/bmag_temp
        bz_normed = bz_temp/bmag_temp

        b_vec = np.array([bx_normed, by_normed, bz_normed])

        #if b_vec.shape == [3,]:
        dot_min = np.dot(b_vec, x_vec)
        if dot_min <= 1 and dot_min >= -1:
            angle = math.degrees (math.acos(dot_min))
            if angle > 90:
                angle = 180 - angle
        else:
            error_list_NaN.append(i)
            angle = np.nan 

        #else:
            #error_list_shape.append(i)
            #angle = np.nan

        omni_cone_angle.append(angle)

    omni_cone = pd.DataFrame({'cone angle': omni_cone_angle, 'datetime': om_df.index})
    omni_cone = omni_cone.set_index('datetime')
    om_df = om_df.join(omni_cone)
    
    print('error no:', len(error_list_NaN), 'total obs:', om_df.shape[0])

    
    return(om_df)

    
