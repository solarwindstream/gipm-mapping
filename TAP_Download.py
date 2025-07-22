import tarfile 
import datetime as dt
#import cdflib
import pandas as pd
from requests import get # to make GET request

#define download function for calling data

def download(url, params, file_name):
    # open in binary mode
    with open(file_name, "wb") as file:
        # get request
        response = get(url, params=params)
        # write to file
        file.write(response.content)
            
#Take interval list and download corresponding tarfiles. Return list of downloaded tarfile filenames.

def dl_Cluster_data_single(int_df, sc_str, filepath):

    #list of interval start and end points
    ints_start_list = int_df['# SC: 1'].tolist()
    index_start = 0
    ints_end_list = int_df[' 2'].tolist()
    
    myurl = 'https://csa.esac.esa.int/csa-sl-tap/data'

    #step one: iterate over every one of the intervals, downloading the tarfiles into a folder
    #note filenames in list for later use
    tarfilelist = []

    if sc_str == 'C1':
        for i, j in zip(ints_start_list, ints_end_list):
            filename = filepath + i + 'C1.tgz'
            tarfilelist.append(filename)
            query_specs = {'RETRIEVAL_TYPE': 'product',
                       'DATASET_ID': 'C1_CP_FGM_FULL',
                       'START_DATE': i,
                       'END_DATE': j,
                       'DELIVERY_FORMAT': 'CDF',
                       'DELIVERY_INTERVAL': 'daily',
                          }
            download(myurl, query_specs, filename)
            
    if sc_str == 'C2':
        for i, j in zip(ints_start_list, ints_end_list):
            filename = filepath + i + 'C2.tgz'
            tarfilelist.append(filename)
            query_specs = {'RETRIEVAL_TYPE': 'product',
                       'DATASET_ID': 'C2_CP_FGM_FULL',
                       'START_DATE': i,
                       'END_DATE': j,
                       'DELIVERY_FORMAT': 'CDF',
                       'DELIVERY_INTERVAL': 'daily'
                          }
            download(myurl, query_specs, filename)
    
    if sc_str == 'C3':
        for i, j in zip(ints_start_list, ints_end_list):
            filename = filepath + i + 'C3.tgz'
            tarfilelist.append(filename)
            query_specs = {'RETRIEVAL_TYPE': 'product',
                       'DATASET_ID': 'C3_CP_FGM_FULL',
                       'START_DATE': i,
                       'END_DATE': j,
                       'DELIVERY_FORMAT': 'CDF',
                       'DELIVERY_INTERVAL': 'daily'
                          }
            download(myurl, query_specs, filename)

    if sc_str == 'C4':
        for i, j in zip(ints_start_list, ints_end_list):
            filename = filepath + i + 'C4.tgz'
            tarfilelist.append(filename)
            query_specs = {'RETRIEVAL_TYPE': 'product',
                       'DATASET_ID': 'C4_CP_FGM_FULL',
                       'START_DATE': i,
                       'END_DATE': j,
                       'DELIVERY_FORMAT': 'CDF',
                       'DELIVERY_INTERVAL': 'daily'
                          }
            download(myurl, query_specs, filename)

    return(tarfilelist)
