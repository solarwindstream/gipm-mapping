import pandas as pd
import datetime as dt

#function to break up Cluster data into 4 min time windows excluding any that do not have full resolution (5280)
#returning a list of start times for each window
#input df is expected to have datetime index

def window_det(df_z):
    
    #set the window size
    time_window = dt.timedelta(seconds=240)
    
    #find datetime of first reading in dataframe
    first_reading = df_z.index[0]

    #(I do not remember why I coded this to be a Series...)
    first_red_s = pd.Series(first_reading)
    #round the first dt reading to the nearest minute
    first_red_s = first_red_s.dt.round('1min')
    first_window = first_red_s[0]
    #also determine the last reading in this dataframe, rounding it to the nearest minute
    length_ind = df_z.index.size
    last_entry = length_ind - 1
    last_reading = df_z.index[last_entry]

    last_red_s = pd.Series(last_reading)
    last_red_s = last_red_s.dt.round('1min')
    last_reading = last_red_s[0]

    #now produce list of possible four minute windows, starting from the rounded start time.
    #...there was really no need to round the end time!

    windows = []
    windows.append(first_window)

    for i in windows:
        if i <= last_reading:
            new_window = i + time_window
            windows.append(new_window)
        else: 
            break
            
    #now, take this list of windows, apply them to the Cluster data, and skip any that do not have data for the whole four minutes
    #full cadence gives 5280 per window. 22*4*60
    #...technically full cadence is 22.4Hz instead of 22Hz, but 22 is high enough resolution for our purposes.
    
    only_full_windows = []

    for i in windows: 
        start_time = i
        end_time = i + time_window

        mask = df_z.loc[(df_z.index >= start_time) & (df_z.index < end_time)]
    
        if mask.shape[0]>=(22*4*60) :
            only_full_windows.append(start_time)
            
    return(only_full_windows)

    
