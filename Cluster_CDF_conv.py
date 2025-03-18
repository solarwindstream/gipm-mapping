#read in CDF file
import cdflib
import pandas as pd
import numpy as np

pd.options.mode.chained_assignment = None

def Cluster_cdf_conv(filename_a, sc_string):
    
    if sc_string=='C1':
        cdf_file_a = cdflib.CDF(filename_a)
        datetimes_a = cdflib.cdfepoch.encode(cdf_file_a['time_tags__C1_CP_FGM_FULL'])

        #convert strings to datetimes
        datetimes_series = pd.Series(datetimes_a) 
        datetimes_a_1 = pd.to_datetime(datetimes_series, errors = 'coerce')
        datetimes_a_1 = datetimes_a_1.to_frame(name="datetime")

        #check that file has data (otherwise causes error)
        shape_a = cdf_file_a['B_vec_xyz_gse__C1_CP_FGM_FULL'].shape

        #shape_a should be length 2. nest whole fnctn w/in this test.
        if len(shape_a) == 2:

            df_a = pd.DataFrame({'Bx_gse': cdf_file_a['B_vec_xyz_gse__C1_CP_FGM_FULL'][:,0], 'By_gse': cdf_file_a['B_vec_xyz_gse__C1_CP_FGM_FULL'][:,1], 
                                   'Bz_gse': cdf_file_a['B_vec_xyz_gse__C1_CP_FGM_FULL'][:,2], "B_mag": cdf_file_a['B_mag__C1_CP_FGM_FULL'], 'X_gse': cdf_file_a['sc_pos_xyz_gse__C1_CP_FGM_FULL'][:,0], 'Y_gse': cdf_file_a['sc_pos_xyz_gse__C1_CP_FGM_FULL'][:,1], 'Z_gse': cdf_file_a['sc_pos_xyz_gse__C1_CP_FGM_FULL'][:,2]})


            df_a = df_a.join(datetimes_a_1)

            #replace fill vals
            df_a[df_a['Bx_gse'] == -1e+30] = np.nan
            df_a[df_a['By_gse'] == -1e+30] = np.nan
            df_a[df_a['Bz_gse'] == -1e+30] = np.nan
            df_a[df_a['B_mag'] == -1e+30] = np.nan

            df_a = df_a.set_index('datetime')


            df_a['R_GSE'] = (df_a['X_gse']**2 + df_a['Y_gse']**2 + df_a['Z_gse']**2)**0.5

            df_a['X_gse'] = df_a['X_gse']/6371
            df_a['Y_gse'] = df_a['Y_gse']/6371
            df_a['Z_gse'] = df_a['Z_gse']/6371
            df_a['R_GSE'] = df_a['R_GSE']/6371

        else:
            print(filename_a, 'data error')
            #return empty df
            df_a = pd.DataFrame({'A' : []})
                
    if sc_string=='C2':
        
        cdf_file_a = cdflib.CDF(filename_a)
        datetimes_a = cdflib.cdfepoch.encode(cdf_file_a['time_tags__C2_CP_FGM_FULL'])

        #convert strings to datetimes
        datetimes_series = pd.Series(datetimes_a) 
        datetimes_a_1 = pd.to_datetime(datetimes_series, errors = 'coerce')
        datetimes_a_1 = datetimes_a_1.to_frame(name="datetime")

        #check that file has data (otherwise causes error)
        shape_a = cdf_file_a['B_vec_xyz_gse__C2_CP_FGM_FULL'].shape

        #shape_a should be length 2. nest whole fnctn w/in this test.
        
        if len(shape_a) == 2:


            df_a = pd.DataFrame({'Bx_gse': cdf_file_a['B_vec_xyz_gse__C2_CP_FGM_FULL'][:,0], 'By_gse': cdf_file_a['B_vec_xyz_gse__C2_CP_FGM_FULL'][:,1], 
                                   'Bz_gse': cdf_file_a['B_vec_xyz_gse__C2_CP_FGM_FULL'][:,2], "B_mag": cdf_file_a['B_mag__C2_CP_FGM_FULL'], 'X_gse': cdf_file_a['sc_pos_xyz_gse__C2_CP_FGM_FULL'][:,0], 'Y_gse': cdf_file_a['sc_pos_xyz_gse__C2_CP_FGM_FULL'][:,1], 'Z_gse': cdf_file_a['sc_pos_xyz_gse__C2_CP_FGM_FULL'][:,2]})


            df_a = df_a.join(datetimes_a_1)

            #replace fill vals
            df_a[df_a['Bx_gse'] == -1e+31] = np.nan
            df_a[df_a['By_gse'] == -1e+31] = np.nan
            df_a[df_a['Bz_gse'] == -1e+31] = np.nan
            df_a[df_a['B_mag'] == -1e+31] = np.nan

            df_a = df_a.set_index('datetime')


            df_a['R_GSE'] = (df_a['X_gse']**2 + df_a['Y_gse']**2 + df_a['Z_gse']**2)**0.5

            df_a['X_gse'] = df_a['X_gse']/6371
            df_a['Y_gse'] = df_a['Y_gse']/6371
            df_a['Z_gse'] = df_a['Z_gse']/6371
            df_a['R_GSE'] = df_a['R_GSE']/6371

        else:
            print(filename_a, 'data error')
            #return empty df
            df_a = pd.DataFrame({'A' : []})
            
    if sc_string=='C3':
        
        cdf_file_a = cdflib.CDF(filename_a)
        datetimes_a = cdflib.cdfepoch.encode(cdf_file_a['time_tags__C3_CP_FGM_FULL'])

        #convert strings to datetimes
        datetimes_series = pd.Series(datetimes_a) 
        datetimes_a_1 = pd.to_datetime(datetimes_series, errors = 'coerce')
        datetimes_a_1 = datetimes_a_1.to_frame(name="datetime")

        #check that file has data (otherwise causes error)
        shape_a = cdf_file_a['B_vec_xyz_gse__C3_CP_FGM_FULL'].shape

        #shape_a should be length 2. nest whole fnctn w/in this test.
        if len(shape_a) == 2:

            df_a = pd.DataFrame({'Bx_gse': cdf_file_a['B_vec_xyz_gse__C3_CP_FGM_FULL'][:,0], 'By_gse': cdf_file_a['B_vec_xyz_gse__C3_CP_FGM_FULL'][:,1], 
                                       'Bz_gse': cdf_file_a['B_vec_xyz_gse__C3_CP_FGM_FULL'][:,2], "B_mag": cdf_file_a['B_mag__C3_CP_FGM_FULL'], 'X_gse': cdf_file_a['sc_pos_xyz_gse__C3_CP_FGM_FULL'][:,0], 'Y_gse': cdf_file_a['sc_pos_xyz_gse__C3_CP_FGM_FULL'][:,1], 'Z_gse': cdf_file_a['sc_pos_xyz_gse__C3_CP_FGM_FULL'][:,2]})


            df_a = df_a.join(datetimes_a_1)

            #replace fill vals
            df_a[df_a['Bx_gse'] == -1e+30] = np.nan
            df_a[df_a['By_gse'] == -1e+30] = np.nan
            df_a[df_a['Bz_gse'] == -1e+30] = np.nan
            df_a[df_a['B_mag'] == -1e+30] = np.nan

            df_a = df_a.set_index('datetime')


            df_a['R_GSE'] = (df_a['X_gse']**2 + df_a['Y_gse']**2 + df_a['Z_gse']**2)**0.5

            df_a['X_gse'] = df_a['X_gse']/6371
            df_a['Y_gse'] = df_a['Y_gse']/6371
            df_a['Z_gse'] = df_a['Z_gse']/6371
            df_a['R_GSE'] = df_a['R_GSE']/6371
    
        else:
            print(filename_a, 'data error')
            #return empty df
            df_a = pd.DataFrame({'A' : []})    
    
    if sc_string=='C4':
        
        cdf_file_a = cdflib.CDF(filename_a)
        datetimes_a = cdflib.cdfepoch.encode(cdf_file_a['time_tags__C4_CP_FGM_FULL'])

        #convert strings to datetimes
        datetimes_series = pd.Series(datetimes_a) 
        datetimes_a_1 = pd.to_datetime(datetimes_series, errors = 'coerce')
        datetimes_a_1 = datetimes_a_1.to_frame(name="datetime")

        #check that file has data (otherwise causes error)
        shape_a = cdf_file_a['B_vec_xyz_gse__C4_CP_FGM_FULL'].shape

        #shape_a should be length 2. nest whole fnctn w/in this test.
        if len(shape_a) == 2:

            df_a = pd.DataFrame({'Bx_gse': cdf_file_a['B_vec_xyz_gse__C4_CP_FGM_FULL'][:,0], 'By_gse': cdf_file_a['B_vec_xyz_gse__C4_CP_FGM_FULL'][:,1], 
                                       'Bz_gse': cdf_file_a['B_vec_xyz_gse__C4_CP_FGM_FULL'][:,2], "B_mag": cdf_file_a['B_mag__C4_CP_FGM_FULL'], 'X_gse': cdf_file_a['sc_pos_xyz_gse__C4_CP_FGM_FULL'][:,0], 'Y_gse': cdf_file_a['sc_pos_xyz_gse__C4_CP_FGM_FULL'][:,1], 'Z_gse': cdf_file_a['sc_pos_xyz_gse__C4_CP_FGM_FULL'][:,2]})


            df_a = df_a.join(datetimes_a_1)

            #replace fill vals
            df_a[df_a['Bx_gse'] == -1e+30] = np.nan
            df_a[df_a['By_gse'] == -1e+30] = np.nan
            df_a[df_a['Bz_gse'] == -1e+30] = np.nan
            df_a[df_a['B_mag'] == -1e+30] = np.nan

            df_a = df_a.set_index('datetime')


            df_a['R_GSE'] = (df_a['X_gse']**2 + df_a['Y_gse']**2 + df_a['Z_gse']**2)**0.5

            df_a['X_gse'] = df_a['X_gse']/6371
            df_a['Y_gse'] = df_a['Y_gse']/6371
            df_a['Z_gse'] = df_a['Z_gse']/6371
            df_a['R_GSE'] = df_a['R_GSE']/6371

        else:
            print(filename_a, 'data error')
            #return empty df
            df_a = pd.DataFrame({'A' : []})

    
    return df_a

    
