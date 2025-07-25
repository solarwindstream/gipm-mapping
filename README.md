# gipm-mapping
Code for processing Cluster and OMNI data to produce GIPM maps of foreshock wave properties

OMNI Data Processing:
1. Download v3.9 CDFs by the year from CDAweb (easier chunks to work with, but not too many) 
2. Run OMNI_cdf_conv.py (makes CSV in correct format to integrate with rest of code) 
3. With list of interval starts, run relevant omni_seg module iteration.

Cluster Data Processing:
1. With input CSVs of relevant datetime intervals from Cluster datamining tool, use Cluster_TAP_Download (in turn using TAP_Download) to obtain tarfiles of relevant Cluster magnetometer and location data from the archive.
2. Run tarfile extract, changing the relevant folder to obtain tgzs from and to export into.

Combined Processing:
1. 


Cluster_CDF_conv: module takes input CDF path and spacecraft number and converts CDF to output dataframe
CDF_CSV: Apocrita script for GIPM conversion.