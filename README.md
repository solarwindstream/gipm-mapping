# gipm-mapping
Code for processing Cluster and OMNI data to produce GIPM maps of foreshock wave properties

OMNI Data Processing:
1. Download CDFs by the year (easier chunks to work with, but not too many) 
2. Run OMNI_cdf_conv.py (makes CSV in correct format to integrate with rest of code) 
3. With list of interval starts, run relevant omni_seg module iteration.
