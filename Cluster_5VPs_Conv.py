#Cluster 5VPs CDF conversion

from Cluster_CDF_conv_5VPs import Cluster_cdf_conv_5VPs

#Cluster_cdf_conv_5VPs(filename_a, sc_string)
scs=['C1', 'C2', 'C3', 'C4']
#filenames=[f'/Users/roseatkinson/Downloads/CSA_Download_20261019_1252/ for sc in scs]
filenames=['/Users/roseatkinson/Downloads/CSA_Download_20261019_1314/C1_CP_FGM_5VPS/C1_CP_FGM_5VPS__20030218_130000_20030218_150000_V140304.cdf', '/Users/roseatkinson/Downloads/CSA_Download_20261019_1314/C2_CP_FGM_5VPS/C2_CP_FGM_5VPS__20030218_130000_20030218_150000_V140305.cdf', '/Users/roseatkinson/Downloads/CSA_Download_20261019_1314/C3_CP_FGM_5VPS/C3_CP_FGM_5VPS__20030218_130000_20030218_150000_V140305.cdf', '/Users/roseatkinson/Downloads/CSA_Download_20261019_1314/C4_CP_FGM_5VPS/C4_CP_FGM_5VPS__20030218_130000_20030218_150000_V140305.cdf']

for fname, sc in zip(filenames,scs):
    df=Cluster_cdf_conv_5VPs(fname, sc)
    firstwin = df.index[0]
    firstwin = str(firstwin)
    CSV_path = '/Users/roseatkinson/Documents/5Hz_Spectra'
    fpath = f"{CSV_path}/{firstwin}_{sc}.csv"
    df.to_csv(fpath)

