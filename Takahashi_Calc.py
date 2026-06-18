#Save Takahashi and Heilig values to Cluster CSV files.

import pandas as pd
import numpy as np
import datetime as dt
import glob

cl_file_list = []

path = "/Users/roseatkinson/Documents/Cluster_Int_CSVs/**"
#test first on duplicated and separate files!
#path = '/Users/roseatkinson/Documents/Testing/**'

for path in glob.glob(path, recursive=True):
    if '.csv' in path:
        cl_file_list.append(path)

for file in cl_file_list:
    df = pd.read_csv(file,encoding='utf-8')
    df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
    df.set_index('datetime', inplace = True)
    
    #Add a row with the Takahashi and Heilig predictions

    #Takahashi prefactor
    pi = np.pi
    m_p = 1.67E-27
    q_p = 1.60E-19
    takahashi_pref = q_p/(4*pi*m_p)

    df['Takahashi Frequency, Hz'] = 1E-9*takahashi_pref*(df['IMF B (mean)'])*(np.cos(np.deg2rad(df['cone angle (mean)'])))**2
    df['Heilig Frequency, Hz'] = 1E-3*(0.78*df['M_A (mean)'] + 0.64)*df['IMF B (mean)']

    df['Takahashi Transverse Difference'] = (df['Takahashi Frequency, Hz'] - df['Peak Transverse Frequency'])
    df['Takahashi Compressive Difference'] = (df['Takahashi Frequency, Hz'] - df['Peak Compressive Frequency'])
    df['Heilig Transverse Difference'] = (df['Heilig Frequency, Hz'] - df['Peak Transverse Frequency'])
    df['Heilig Compressive Difference'] = (df['Heilig Frequency, Hz'] - df['Peak Compressive Frequency'])

    df.to_csv(file)

