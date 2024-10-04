import tarfile 
import datetime as dt
#import cdflib
import pandas as pd
from requests import get # to make GET request
import numpy as np
import glob

#define download function for calling data
batch_no = 1

def download(url, params, file_name):
    # open in binary mode
    with open(file_name, "wb") as file:
        # get request
        response = get(url, params=params)
        # write to file
        file.write(response.content)
            
#input data list and download tarfiles. return list of filenames

def dl_Cluster_data_single(int_df, sc_no):

    #list of interval start and end points
    ints_start_list = int_df['# SC: 1'].tolist()
    index_start = 0
    ints_end_list = int_df[' 2'].tolist()
    
    myurl = 'https://csa.esac.esa.int/csa-sl-tap/data'

    #step one: iterate over every one of the intervals, downloading the tarfiles into a folder
    #note filenames in list for later use
    tarfilelist = []

    if sc_no == '1':
        for i, j in zip(ints_start_list, ints_end_list):
            filename = '/data/scratch/apx059/23_Years_Data/C1/' + i + '.tgz'
            tarfilelist.append(filename)
            query_specs = {'RETRIEVAL_TYPE': 'product',
                       'DATASET_ID': 'C1_CP_FGM_FULL',
                       'START_DATE': i,
                       'END_DATE': j,
                       'DELIVERY_FORMAT': 'CDF',
                       'DELIVERY_INTERVAL': 'daily',
                          }
            download(myurl, query_specs, filename)
            
    if sc_no == '2':
        for i, j in zip(ints_start_list, ints_end_list):
            filename = '/data/scratch/apx059/23_Years_Data/C2/' + i + '.tgz'
            tarfilelist.append(filename)
            query_specs = {'RETRIEVAL_TYPE': 'product',
                       'DATASET_ID': 'C2_CP_FGM_FULL',
                       'START_DATE': i,
                       'END_DATE': j,
                       'DELIVERY_FORMAT': 'CDF',
                       'DELIVERY_INTERVAL': 'daily'}
            download(myurl, query_specs, filename)
    
    if sc_no == '3':
        for i, j in zip(ints_start_list, ints_end_list):
            filename = '/data/scratch/apx059/23_Years_Data/C3/' + i + '.tgz'
            tarfilelist.append(filename)
            query_specs = {'RETRIEVAL_TYPE': 'product',
                       'DATASET_ID': 'C3_CP_FGM_FULL',
                       'START_DATE': i,
                       'END_DATE': j,
                       'DELIVERY_FORMAT': 'CDF',
                       'DELIVERY_INTERVAL': 'daily'}
            download(myurl, query_specs, filename)

    if sc_no == '4':
        for i, j in zip(ints_start_list, ints_end_list):
            filename = '/data/scratch/apx059/23_Years_Data/C4/' + i + '.tgz'
            tarfilelist.append(filename)
            query_specs = {'RETRIEVAL_TYPE': 'product',
                       'DATASET_ID': 'C4_CP_FGM_FULL',
                       'START_DATE': i,
                       'END_DATE': j,
                       'DELIVERY_FORMAT': 'CDF',
                       'DELIVERY_INTERVAL': 'daily'}
            download(myurl, query_specs, filename)

    return(tarfilelist)

csv_name = '/data/scratch/apx059/Cluster_Intervals-01022002-01022024/dm-intervals-c2-240902-144411.csv'

cluster_2_23yr = pd.read_csv(csv_name)

#drop unnecessary rows (without data in!)
un_rows = [0,1,2,3,4]
cluster_2_23yr = cluster_2_23yr.drop(un_rows)
cluster_2_23yr = cluster_2_23yr.reindex()

tf_list = dl_Cluster_data_single(cluster_2_23yr, '2')
