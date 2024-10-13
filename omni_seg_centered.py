#averaging time window of ten minutes, centered around midpoint of Cluster interval (bearing in mind ref time is window start)
#this function takes 1 min OMNI readings and a list of time intervals and produces df of omni average B, V and n values
#including components
#for 10 minute window centered on midpoint of time interval
import datetime as dt
import pandas as pd
import numpy as np
import math


def omni_seg(om_df, only_full_windows):
    time_back = dt.timedelta(minutes=4)
    time_forward = dt.timedelta(minutes=6)

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
        
        
    om_averages = pd.DataFrame({'datetime': only_full_windows, 'Np': omni_N_ave, 'B_mag': omni_B_ave, 'V_gse': omni_V_ave, 'B_X_gse': omni_Bx_ave, 'B_Y_gse': omni_By_ave, 'B_Z_gse': omni_Bz_ave, 'V_X_gse': omni_Vx_ave, 'V_Y_gse': omni_Vy_ave, 'V_Z_gse': omni_Vz_ave, 'M_A':omni_MA_ave, 'cone angle':omni_CA_ave})

    om_averages['Ave B']= (om_averages['B_X_gse']**2 + om_averages['B_Y_gse']**2 + om_averages['B_Z_gse']**2)**0.5
    om_averages['Norm Bx'] = om_averages['B_X_gse']/om_averages['Ave B']
    om_averages['Norm By'] = om_averages['B_Y_gse']/om_averages['Ave B']
    om_averages['Norm Bz'] = om_averages['B_Z_gse']/om_averages['Ave B']
    
    return(om_averages)

