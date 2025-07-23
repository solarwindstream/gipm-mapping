#Run TAP module

from TAP_Download import download, dl_Cluster_data_single
import pandas as pd

csv_1 = '/data/home/apx059/Cluster_Intervals_01022002-30092024/dm-intervals-c1-250722-161131.csv'
csv_2 = '/data/home/apx059/Cluster_Intervals_01022002-30092024/dm-intervals-c2-250722-161131.csv'
csv_3 = '/data/home/apx059/Cluster_Intervals_01022002-30092024/dm-intervals-c3-250722-161131.csv'
csv_4 = '/data/home/apx059/Cluster_Intervals_01022002-30092024/dm-intervals-c4-250722-161131.csv'

cluster_1_23 = pd.read_csv(csv_1)
cluster_2_23 = pd.read_csv(csv_2)
cluster_3_23 = pd.read_csv(csv_3)
cluster_4_23 = pd.read_csv(csv_4)

#drop unnecessary rows (without data in!)
un_rows = [0,1,2,3,4]

cluster_1_23 = cluster_1_23.drop(un_rows)
cluster_1_23 = cluster_1_23.reindex()

cluster_2_23 = cluster_2_23.drop(un_rows)
cluster_2_23 = cluster_2_23.reindex()

cluster_3_23 = cluster_3_23.drop(un_rows)
cluster_3_23 = cluster_3_23.reindex()

cluster_4_23 = cluster_4_23.drop(un_rows)
cluster_4_23 = cluster_4_23.reindex()

tf_list_1 = dl_Cluster_data_single(cluster_1_23, 'C1', '/data/SPCS-HIETALA-Shocks/GIPM-MAPPING/Cluster_Source_Tarfiles/')

tf_list_2 = dl_Cluster_data_single(cluster_2_23, 'C2', '/data/SPCS-HIETALA-Shocks/GIPM-MAPPING/Cluster_Source_Tarfiles/')

tf_list_3 = dl_Cluster_data_single(cluster_3_23, 'C3', '/data/SPCS-HIETALA-Shocks/GIPM-MAPPING/Cluster_Source_Tarfiles/')

tf_list_4 = dl_Cluster_data_single(cluster_4_23, 'C4', '/data/SPCS-HIETALA-Shocks/GIPM-MAPPING/Cluster_Source_Tarfiles/')

#dl_Cluster_data_single(int_df, sc_str, filepath) returns tarfilelist

#will be saving tarfiles into /data/SPCS-Hietala-Shocks/GIPM-MAPPING/Cluster_Source_Tarfiles