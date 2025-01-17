#make histogram plots for all 23 years

import pandas as pd
import numpy as np
import datetime as dt
import glob

from XMA_finder import XMA_finder
from histo_plot import histo_plot
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

#import modules
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib import cm
import matplotlib
from merka05_surface_eq_array_GIPM import merka05_surface_eq_array_GIPM

#open OMNI CSVs for all 23 years. NOT the averaged ones. 

omni_path = '/data/scratch/apx059/OMNI_Raw/' 
omni_csvs=['omni_hros_1min_20010201000000_20020201000000.csv','omni_hros_1min_20020201000000_20030201000000.csv','omni_hros_1min_20030201000000_20040201000000.csv','omni_hros_1min_20040201000000_20050201000000.csv','omni_hros_1min_20050201000000_20060201000000.csv','omni_hros_1min_20060201000000_20070201000000.csv','omni_hros_1min_20070201000000_20080201000000.csv','omni_hros_1min_20080201000000_20090201000000.csv','omni_hros_1min_20090201000000_20100201000000.csv','omni_hros_1min_20100201000000_20110201000000.csv','omni_hros_1min_20110201000000_20120201000000.csv','omni_hros_1min_20120201000000_20130201000000.csv','omni_hros_1min_20130201000000_20140201000000.csv','omni_hros_1min_20140201000000_20150201000000.csv','omni_hros_1min_20150201000000_20160201000000.csv','omni_hros_1min_20160201000000_20170201000000.csv','omni_hros_1min_20170201000000_20180201000000.csv','omni_hros_1min_20180201000000_20190201000000.csv','omni_hros_1min_20190201000000_20200201000000.csv','omni_hros_1min_20200201000000_20210201000000.csv','omni_hros_1min_20210201000000_20220201000000.csv','omni_hros_1min_20220201000000_20230201000000.csv','omni_hros_1min_20230201000000_20240201000000.csv']

for element in omni_csvs:
    if '.csv' in element:
        filepath = omni_path + element
        om_csvs.append(filepath)
        
om_dfs = []

for element in om_csvs:
    om = pd.read_csv(element)
    om['datetime'] = pd.to_datetime(om['datetime'])
    om = om.set_index('datetime', inplace = True)
    om_dfs.append(om)
    
omni_all = pd.concat(om_dfs)

XMA_all = XMA_finder(omni_all)

##shouldn't the XMA really be the one used in the Merka model?

