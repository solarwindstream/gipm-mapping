import pandas as pd
import numpy as np
import datetime as dt
import glob

from Cluster_CDF_conv import Cluster_cdf_conv

##load Cluster CDFs

cl_file_list_C1 = []
cl_file_list_C2 = []
cl_file_list_C3 = []
cl_file_list_C4 = []

path = "/Users/roseatkinson/Documents/Cluster_Time_Series/**"

for path in glob.glob(path, recursive=True):
    if '.cdf' in path:
        if 'C1' in path:
            cl_file_list_C1.append(path)
        if 'C2' in path:
            cl_file_list_C2.append(path)
        if 'C3' in path:
            cl_file_list_C3.append(path)
        if 'C4' in path:
            cl_file_list_C4.append(path)

    
#load in cdfs & create list of dataframes of raw Cluster data from CDF file list!
    
err_list: list[str] = []

#convert to CSV and save as CSV.

for count, file in enumerate(cl_file_list_C1):
    df = Cluster_cdf_conv(file, 'C1')
    if not df.empty:
        first_time = df.index[0]
        time_str = first_time.strftime("%Y-%m-%d %H:%M:%S")
        csv_name = 'TS_'+time_str+'_C1.csv'
        csv_path = "/Users/roseatkinson/Documents/Cluster_Time_Series/"+csv_name
        df.to_csv(csv_path)
    else:
        err_list.append(file)
    print(count+1, 'out of', len(cl_file_list_C1), 'C1 CDFs loaded in', flush=True)


for count, file in enumerate(cl_file_list_C2):
    df = Cluster_cdf_conv(file, 'C2')
    if not df.empty:
        first_time = df.index[0]
        time_str = first_time.strftime("%Y-%m-%d %H:%M:%S")
        csv_name = 'TS_'+time_str+'_C2.csv'
        csv_path = "/Users/roseatkinson/Documents/Cluster_Time_Series/"+csv_name
        df.to_csv(csv_path)
    else:
        err_list.append(file)
    print(count+1, 'out of', len(cl_file_list_C2), 'C2 CDFs loaded in', flush=True)
    
for count, file in enumerate(cl_file_list_C3):
    df = Cluster_cdf_conv(file, 'C3')
    if not df.empty:
        first_time = df.index[0]
        time_str = first_time.strftime("%Y-%m-%d %H:%M:%S")
        csv_name = 'TS_'+time_str+'_C3.csv'
        csv_path = "/Users/roseatkinson/Documents/Cluster_Time_Series/"+csv_name
        df.to_csv(csv_path)
    else:
        err_list.append(file)
    print(count+1, 'out of', len(cl_file_list_C3), 'C3 CDFs loaded in', flush=True)

for count, file in enumerate(cl_file_list_C4):
    df = Cluster_cdf_conv(file, 'C4')
    if not df.empty:
        first_time = df.index[0]
        time_str = first_time.strftime("%Y-%m-%d %H:%M:%S")
        csv_name = 'TS_'+time_str+'_C4.csv'
        csv_path = "/Users/roseatkinson/Documents/Cluster_Time_Series/"+csv_name
        df.to_csv(csv_path)
    else:
        err_list.append(file)
    print(count+1, 'out of', len(cl_file_list_C4), 'C4 CDFs loaded in', flush=True)