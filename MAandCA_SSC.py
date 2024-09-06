#Split location-only results by MA and cone angle, and save as new CSVs

import pandas as pd
import numpy as np
import datetime as dt
import glob

batch = 0
batch_no = 'batch' + str(batch)
start_file = batch*5
end_file = (batch+1)*5

#open OMNI *average* CSVs
list_all = []

path = "/Users/apx059/Documents/Location_Only_Checks/OMNI Averages/**"

for path in glob.glob(path, recursive=True):
    list_all.append(path)
    
#list with only files, not folders
om_csvs = []

for element in list_all:
    if '.csv' in element:
        om_csvs.append(element)
        
om_dfs = []

for element in om_csvs:
    om = pd.read_csv(element)
    om_dfs.append(om)
    
omni_all = pd.concat(om_dfs)
omni_all['datetime'] = pd.to_datetime(omni_all['datetime'])

omni_all = omni_all.set_index('datetime')

for element in om_dfs:
    element['datetime'] = pd.to_datetime(element['datetime'])
    element = element.set_index('datetime', inplace = True)

#load Cluster CSVs

list_all = []

path = "/Users/apx059/Documents/OM_SEG_TEST/temp_cl_dfs/**"

for path in glob.glob(path, recursive=True):
    list_all.append(path)
    
#list with only files, not folders
cl_file_list = []

for element in list_all:
    if '.csv' in element:
        cl_file_list.append(element)

cl_dfs = []

for file in cl_file_list[0:5]:
    df = pd.read_csv(file,encoding='utf-8', names = ['datetime', 'GIPM X (RE)', 'GIPM Y (RE)', 'GIPM Z (RE)'])
    df['datetime'] = pd.to_datetime(df['datetime'])
    df.set_index('datetime', inplace = True)
    cl_dfs.append(df)


#iterate over all cl_dfs 

for df in cl_dfs:
    
    ma_list = []
    cone_a_list = []
    
    for i in df.index:
        if i in omni_all.index:
            ma_temp = omni_all['M_A'].loc[i]
            cone_a_temp = omni_all['cone angle'].loc[i]
            ma_list.append(ma_temp)
            cone_a_list.append(cone_a_temp)
        else:
            ma_list.append(np.nan)
            cone_a_list.append(np.nan)
            
    df['MA'] = ma_list
    df['cone angle'] = cone_a_list

#save results as CSVs

csvpath = "/Users/apx059/Documents/OM_SEG_TEST/"

startref = 1

for df in cl_dfs:
    #saving file
    ref = str(startref)
    start_year = df.index[0].strftime("%Y")
    start_year_int = int(start_year)
    end_year_int = start_year_int + 1
    end_year = str(end_year_int)
    filename = csvpath + 'Feb' + start_year + '_Feb' + end_year + 'df_' + ref + batch_no + '.csv'
    df.to_csv(filename)
    startref = startref + 1