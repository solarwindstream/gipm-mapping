#Script to do linear correlation of two time series

from scipy import stats
import datetime as dt
import pandas as pd

import matplotlib 
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.gridspec import GridSpec

from gipm_plotting_modules import model_calcs, draw_background, set_limits

############################
############################

def time_series_comp(df_dict:dict, sc_list:list[str],fname:str) -> None:
    """Take time series dict and spacecraft reference to plot two ts side-by-side"""
    if len(sc_list) > 2:
        print("Spacecraft list too long, only handling first two")
        sc_list=sc_list[:2]
    
    fig = plt.figure(figsize=(8, 3))
        
    locator = mdates.AutoDateLocator(minticks=3, maxticks=5)
    formatter = mdates.ConciseDateFormatter(locator)
        
    # Define a gridspec layout
    gs = GridSpec(1, 2, figure=fig)
    
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    
    for ax, sc in zip([ax1,ax2], sc_list):
        df=df_dict[sc]
        ax.plot(df.index, df['Bx_gse'], color = 'red', label='Bx')
        ax.plot(df.index, df['By_gse'], color = 'green', label='By')
        ax.plot(df.index, df['Bz_gse'], color = 'blue', label='Bz')
        ax.plot(df.index, df['B_mag'], color = 'k', label='Bmag')
        ax.set_ylim(-10, 10)
        ax.set_xlim(df.index[0], df.index[-1])
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.xaxis.set_tick_params(labelsize=12)
        ax.yaxis.set_tick_params(labelsize=12)
        ax.set_ylabel(r'B, nT')
        ax.set_title(sc)
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='lightgray')
        ax.xaxis.grid(color='lightgray')
        if sc==sc_list[1]:
            ax.legend(loc='upper right')
    
    path = f"/Users/roseatkinson/Documents/New_Figs/Correlation_Study/{fname}.png"
            
    plt.savefig(path, bbox_inches='tight')

def loc_plot(locs_x:list[float], locs_y:list[float], fname:str) -> None: 
    """Plots two locations in X,Y GIPM"""
    fig = plt.figure(figsize=(4, 4))
    ax=fig.add_axes((0, 0, 1, 1))
    
    draw_background(ax, Xgipm[:, :, 0], Ygipm[:, :, 0], f[:, :, 0],
                        X_shue, R_shue)
    
    ax.scatter(locs_x, locs_y, marker="o", c=['k','indigo'], s=60)
    
    set_limits(ax)
    
    ax.set_title(r'Case Study Locations')
    ax.set_xlabel("$X_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")
    ax.set_ylabel("$Y_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")
    
    path = f"/Users/roseatkinson/Documents/New_Figs/Correlation_Study/{fname}.png"
    plt.savefig(path, bbox_inches='tight')


############################
############################

int_csv_files =['/Users/roseatkinson/Documents/Cluster_Int_CSVs/2017-03-07 21:18:00_C1.csv','/Users/roseatkinson/Documents/Cluster_Int_CSVs/2017-03-07 22:34:00_C3.csv', '/Users/roseatkinson/Documents/Cluster_Int_CSVs/2017-03-07 22:06:00_C4.csv']

ts_files = ['/Users/roseatkinson/Documents/Cluster_Time_Series/TS_2017-03-07 21:00:00_C1.csv','/Users/roseatkinson/Documents/Cluster_Time_Series/TS_2017-03-07 21:00:00_C3.csv','/Users/roseatkinson/Documents/Cluster_Time_Series/TS_2017-03-07 21:00:00_C4.csv']

df_dict = {}

scs=['C1','C3','C4']

for sc, int_file, ts_file in zip(scs, int_csv_files,ts_files):
    df_dict[sc]={}
    df = pd.read_csv(int_file,encoding='utf-8')
    df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
    df.set_index('datetime', inplace = True)
    df['sc']=sc
    df_dict[sc]['Interval Data']=df
    
    df = pd.read_csv(ts_file,encoding='utf-8')
    df['datetime'] = pd.to_datetime(df['datetime'], format='mixed')
    df.set_index('datetime', inplace = True)
    df_dict[sc]['Time Series']=df

#int_start = datetime.strftime("%Y-%m-%d %H:%M:%S")

int_start = pd.to_datetime("2017-03-08 00:50:00")
int_end = int_start + dt.timedelta(minutes=4)

#for now ignore the integrated data and just find the three matching time series

ts_short={}

for sc in scs:
    ts_df=df_dict[sc]['Time Series']
    ts_short[sc] = ts_df.loc[((ts_df.index >= int_start) & (ts_df.index < int_end))]
    
#Bx_gse, By_gse, Bz_gse, B_mag
for sc in scs: print (ts_short[sc].index[0])

C1_mag = ts_short['C1']['B_mag'].to_numpy()
C3_mag = ts_short['C3']['B_mag'].to_numpy()
C4_mag = ts_short['C4']['B_mag'].to_numpy()

if C3_mag.shape != C1_mag.shape:
    #print('C1, C3 len', C1_mag.shape, C3_mag.shape)
    length = C1_mag.shape[0]
    C3_mag = C3_mag[:length]

if C4_mag.shape != C1_mag.shape:
    #print('C1, C4 len', C1_mag.shape, C4_mag.shape)
    length = C1_mag.shape[0]
    C4_mag = C4_mag[:length]

r_1_1, p_1_1=stats.pearsonr(C1_mag, C1_mag)
r_1_3, p_1_3=stats.pearsonr(C1_mag, C3_mag)
r_1_4, p_1_4=stats.pearsonr(C1_mag, C4_mag)

print(f"test correlation is {r_1_1}, C1-C3 correlation is {r_1_3}, C1-C4 correlation is {r_1_4}")

#get Vsw data, GSE loc for 00:50 C1 case
#'GSE X (mean)' 'SW Vx (mean)'

int_locs = {}

for sc in scs:
    int_locs[sc]={}
    int_locs[sc]['Data']=df_dict[sc]['Interval Data'].loc[[int_start]]
    print(sc)
    for axis in ['GSE X (mean)', 'GSE Y (mean)', 'GSE Z (mean)']:
        print(f"{axis}: {int_locs[sc]['Data'].at[int_start, axis]}")
            
#print(C1_0050_data, C3_0050_data, C4_0050_data)
# print(C1_0050_data.at[int_start, 'GSE X (mean)'], C3_0050_data.at[int_start, 'GSE X (mean)'], C4_0050_data.at[int_start, 'GSE X (mean)'])
#print(C3_0050_data.loc[['GSE X (mean)', 'GSE Y (mean)', 'GSE Z (mean)']])
#print(C4_0050_data.loc[['GSE X (mean)', 'GSE Y (mean)', 'GSE Z (mean)']])

#get GSE loc for 00:50 C3, C4 cases 
#determine if going backward in time or forward in time from C1 by looking (mostly) at GSE X values.
#calculate expected advection time and then shift in a +-30s window around it to find max correlation value

# Create a figure for two instances overlaid
fig = plt.figure(figsize=(5, 3))
ax=fig.add_axes((0, 0, 1, 1))

locator = mdates.AutoDateLocator(minticks=3, maxticks=5)
formatter = mdates.ConciseDateFormatter(locator)

df1=ts_short['C1']
df3=ts_short['C3']

Cmin1 = (df1['B_mag'].min()) - 1
Cmax1 = (df1['B_mag'].max()) + 1

ax.plot(df1.index, df1['B_mag'], color = 'k', label='Bmag_C1')
ax.plot(df3.index, df3['B_mag'], color = 'indigo', label='Bmag_C3')
ax.set_ylim(Cmin1, Cmax1)
ax.set_xlim(df.index[0], df.index[-1])
ax.xaxis.set_major_formatter(formatter)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)
ax.legend(loc='lower right')

path = "/Users/roseatkinson/Documents/New_Figs/Correlation_Study/C1C3overlaid.png"
        
plt.savefig(path, bbox_inches='tight')

X_shue, R_shue, Xgipm, Ygipm, Zgipm, f, fitting_coeffs = model_calcs()

loc_x = []
loc_y = []
loc_z = []

for sc in scs:
    loc_x.append(int_locs[sc]['Data'].at[int_start, 'GIPM X (OMNI mean)'])
    loc_y.append(int_locs[sc]['Data'].at[int_start, 'GIPM Y (OMNI mean)'])
    loc_z.append(int_locs[sc]['Data'].at[int_start, 'GIPM Z (OMNI mean)'])

time_series_comp()
loc_plot(loc_x[0:2],loc_y[0:2],fname='C1C3_2017_03_08_00_50')   

###################

#don't actually need fs!
#fs['/Users/roseatkinson/Documents/Extreme_Spectra/FS_2017-03-08 00:50:00_C1.csv',
# '/Users/roseatkinson/Documents/Extreme_Spectra/FS_2017-03-08 00:50:00_C3.csv',
# '/Users/roseatkinson/Documents/Extreme_Spectra/FS_2017-03-08 00:50:00_C4.csv']
