#split cluster data into cone angle times
#now make dfs for only those date times

#average Cluster readings over each two minute interval
####### this is bit I need to combine for mag readings.

#need to have associated the windows with the OMNI averages

def cone_angle_dfs(om_averages, cluster_df, cluster_dt_loc):
    import datetime as dt
    import pandas as pd
    
    winds_030 = []
    winds_3060 = []
    winds_6090 = []
    
    for i in om_averages.index:
        
        if om_averages.loc[i, 'cone angle'] <= 30:
            winds_030.append(i)
            
        elif om_averages.loc[i, 'cone angle'] <= 60:
            winds_3060.append(i)
            
        else:
            winds_6090.append(i)
    
    time_window = dt.timedelta(seconds=120)
    B_list = om_averages['B_mag'].tolist()
    
    cluster_locs_030 = []
    cluster_locs_3060 = []
    cluster_locs_6090 = []
    cluster_b_rat_min_030 = []
    cluster_b_rat_max_030 = []
    cluster_b_rat_mean_030 = []
    cluster_b_rat_min_3060 = []
    cluster_b_rat_max_3060 = []
    cluster_b_rat_mean_3060 = []
    cluster_b_rat_min_6090 = []
    cluster_b_rat_max_6090 = []
    cluster_b_rat_mean_6090 = []
    
    for i in winds_030:
        start_time = i
        end_time = i + time_window
        mask = cluster_df.loc[(cluster_df.index >= start_time) & (cluster_df.index < end_time)]
        Cluster_list = mask['B_mag'].tolist()
        Cluster_min = min(Cluster_list)
        Cluster_max = max(Cluster_list)
        Cluster_mean = sum(Cluster_list)/len(Cluster_list)
        Omni_ave = om_averages.loc[i,'B_mag']
        ratio_min = Cluster_min/Omni_ave
        ratio_max = Cluster_max/Omni_ave
        ratio_mean = Cluster_mean/Omni_ave
        cluster_b_rat_min_030.append(ratio_min)
        cluster_b_rat_max_030.append(ratio_max)
        cluster_b_rat_mean_030.append(ratio_mean)
        
        j = cluster_dt_loc.index[cluster_dt_loc['datetime']==i]
        cluster_locs_030.append(cluster_dt_loc.loc[j[0],'GIPM Loc'])
    
    for i in winds_3060:
        start_time = i
        end_time = i + time_window
        mask = cluster_df.loc[(cluster_df.index >= start_time) & (cluster_df.index < end_time)]
        Cluster_list = mask['B_mag'].tolist()
        Cluster_min = min(Cluster_list)
        Cluster_max = max(Cluster_list)
        Cluster_mean = sum(Cluster_list)/len(Cluster_list)
        Omni_ave = om_averages.loc[i,'B_mag']
        ratio_min = Cluster_min/Omni_ave
        ratio_max = Cluster_max/Omni_ave
        ratio_mean = Cluster_mean/Omni_ave
        cluster_b_rat_min_3060.append(ratio_min)
        cluster_b_rat_max_3060.append(ratio_max)
        cluster_b_rat_mean_3060.append(ratio_mean)
        
        j = cluster_dt_loc.index[cluster_dt_loc['datetime']==i]
        cluster_locs_3060.append(cluster_dt_loc.loc[j[0],'GIPM Loc'])
        
    for i in winds_6090:
        start_time = i
        end_time = i + time_window
        mask = cluster_df.loc[(cluster_df.index >= start_time) & (cluster_df.index < end_time)]
        Cluster_list = mask['B_mag'].tolist()
        Cluster_min = min(Cluster_list)
        Cluster_max = max(Cluster_list)
        Cluster_mean = sum(Cluster_list)/len(Cluster_list)
        Omni_ave = om_averages.loc[i,'B_mag']
        ratio_min = Cluster_min/Omni_ave
        ratio_max = Cluster_max/Omni_ave
        ratio_mean = Cluster_mean/Omni_ave
        cluster_b_rat_min_6090.append(ratio_min)
        cluster_b_rat_max_6090.append(ratio_max)
        cluster_b_rat_mean_6090.append(ratio_mean)
        
        j = cluster_dt_loc.index[cluster_dt_loc['datetime']==i]
        cluster_locs_6090.append(cluster_dt_loc.loc[j[0],'GIPM Loc'])
        
    x_list_rad = []
    y_list_rad = []
    z_list_rad = []

    x_list_spir = []
    y_list_spir = []
    z_list_spir = []

    x_list_obl = []
    y_list_obl = []
    z_list_obl = []
    
    length = len(cluster_locs_030)
    for i in range(length):
        x_list_rad.append(cluster_locs_030[i][0])
    

    for i in range(length):
        y_list_rad.append(cluster_locs_030[i][1])
    

    for i in range(length):
        z_list_rad.append(cluster_locs_030[i][2])
    
    length = len(cluster_locs_3060)
    for i in range(length):
        x_list_spir.append(cluster_locs_3060[i][0])
    

    for i in range(length):
        y_list_spir.append(cluster_locs_3060[i][1])
    

    for i in range(length):
        z_list_spir.append(cluster_locs_3060[i][2])
    
    length = len(cluster_locs_6090)
    for i in range(length):
        x_list_obl.append(cluster_locs_6090[i][0])
    

    for i in range(length):
        y_list_obl.append(cluster_locs_6090[i][1])
    

    for i in range(length):
        z_list_obl.append(cluster_locs_6090[i][2])

        
    dict_030 = {'window start': winds_030, 'Cluster Loc GIPM X': x_list_rad, 'Cluster Loc GIPM Y': y_list_rad,'Cluster Loc GIPM Z': z_list_rad, 'Bc_Bo Min': cluster_b_rat_min_030, 'Bc_Bo Mean': cluster_b_rat_mean_030, 'Bc_Bo Max':cluster_b_rat_max_030}
    dict_3060 = {'window start': winds_3060, 'Cluster Loc GIPM X': x_list_spir, 'Cluster Loc GIPM Y': y_list_spir,'Cluster Loc GIPM Z': z_list_spir, 'Bc_Bo Min': cluster_b_rat_min_3060, 'Bc_Bo Mean': cluster_b_rat_mean_3060, 'Bc_Bo Max':cluster_b_rat_max_3060}
    dict_6090 = {'window start': winds_6090, 'Cluster Loc GIPM X': x_list_obl, 'Cluster Loc GIPM Y': y_list_obl,'Cluster Loc GIPM Z': z_list_obl, 'Bc_Bo Min': cluster_b_rat_min_6090, 'Bc_Bo Mean': cluster_b_rat_mean_6090, 'Bc_Bo Max':cluster_b_rat_max_6090}
    
    df_030 = pd.DataFrame(dict_030)
    df_3060 = pd.DataFrame(dict_3060)
    df_6090 = pd.DataFrame(dict_6090)   
    
    return (df_030, df_3060, df_6090)
