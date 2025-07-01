#take all upstream conditions and bin into B-V hist
#probably with cone angle split too?

import pandas as pd
import numpy as np
import datetime as dt
import glob

from XMA_finder import XMA_finder
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

#import modules
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib
from merka05_surface_eq_array_GIPM import merka05_surface_eq_array_GIPM
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

#matplotlib.pyplot.hist2d(x, y, bins=10, *, range=None, density=False, weights=None, cmin=None, cmax=None, data=None, **kwargs)

#import *OMNI* data!

omni_all_nodup = omni_all.drop_duplicates()

x = omni_all_nodup['B_mag']
y = omni_all_nodup['V_gse']

fig, ax = plt.subplots()

h = ax.hist2d(x, y, bins=40, range=[[0,20],[200,800]], cmap='Blues', cmin=5000)
ax.set_xlabel(r'$B$, nT')
ax.set_ylabel(r'$V_{sw}$, km/s')
fig.colorbar(h[3], ax=ax)
plt.show()
