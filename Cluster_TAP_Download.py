#Run TAP module

from TAP_Download import download, dl_Cluster_data_single
import pandas as pd

csv_1 = '/data/home/apx059/Cluster_Intervals_16032001-01022002/dm-intervals-c1-240731-131130.csv'

cluster_2001 = pd.read_csv(csv_name)

#drop unnecessary rows (without data in!)
un_rows = [0,1,2,3,4]
cluster_2001 = cluster_2001.drop(un_rows)
cluster_2001 = cluster_2001.reindex()

tf_list = dl_Cluster_data_single(cluster_2001, 'C1', '/data/SPCS-Hietala-Shocks/GIPM-MAPPING/Cluster_Source_Tarfiles/')

#dl_Cluster_data_single(int_df, sc_str, filepath) returns tarfilelist

#will be saving tarfiles into /data/SPCS-Hietala-Shocks/GIPM-MAPPING/Cluster_Source_Tarfiles