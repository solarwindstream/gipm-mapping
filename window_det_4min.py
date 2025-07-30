#function to break up Cluster data into 4 min time windows
#excluding any that do not have full resolution (5280)

def window_det(df_z):
    
    import pandas as pd
    import datetime as dt
    
    #set the window size
    time_window = dt.timedelta(seconds=240)

    first_reading = df_z.index[0]

    first_red_s = pd.Series(first_reading)

    first_red_s = first_red_s.dt.round('1min')
    first_window = first_red_s[0]

    length_ind = df_z.index.size
    last_entry = length_ind - 1
    last_reading = df_z.index[last_entry]

    last_red_s = pd.Series(last_reading)
    last_red_s = last_red_s.dt.round('1min')
    last_reading = last_red_s[0]

    windows = []
    windows.append(first_window)

    for i in windows:
        if i <= last_reading:
            new_window = i + time_window
            windows.append(new_window)
        else: 
            break
            
    #skip any that do not have data for the whole four minutes
    #full cadence gives 5280 per window. 22*4*60
    
    only_full_windows = []

    for i in windows: 
        start_time = i
        end_time = i + time_window

        mask = df_z.loc[(df_z.index >= start_time) & (df_z.index < end_time)]
    
        if mask.shape[0]>=(22*4*60) :
            only_full_windows.append(start_time)
            
    return(only_full_windows)

    
