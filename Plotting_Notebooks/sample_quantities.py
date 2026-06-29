import pandas as pd

def sample_filter(df_high, df_low):

    #filter dfs for FS, SW, QPARA MSH & QPERP MSH
    #for spiral cone angles, 30-45 degrees
    
    #Foreshock
    
    x_lim_fs = [10, 12]
    y_lim_fs = [-13, -11]
    
    #Solar Wind
    
    x_lim_sw = [10, 12]
    y_lim_sw = [11, 13]
    
    
    #Qpara Magnetosheath
    
    x_lim_para_msh = [6, 8]
    y_lim_para_msh = [-13, -11]
    
    
    #Qperp Magnetosheath
    
    x_lim_perp_msh = [6, 8]
    y_lim_perp_msh = [11, 13]

    
    fs_high = df_high.loc[(df_high['GIPM X (OMNI mean)'] > x_lim_fs[0]) & (df_high['GIPM X (OMNI mean)'] < x_lim_fs[1]) & (df_high['GIPM Y (OMNI mean)'] > y_lim_fs[0]) & (df_high['GIPM Y (OMNI mean)'] < y_lim_fs[1])]
    fs_low = df_low.loc[(df_low['GIPM X (OMNI mean)'] > x_lim_fs[0]) & (df_low['GIPM X (OMNI mean)'] < x_lim_fs[1]) & (df_low['GIPM Y (OMNI mean)'] > y_lim_fs[0]) & (df_low['GIPM Y (OMNI mean)'] < y_lim_fs[1])]
    
    sw_high = df_high.loc[(df_high['GIPM X (OMNI mean)'] > x_lim_sw[0]) & (df_high['GIPM X (OMNI mean)'] < x_lim_sw[1]) & (df_high['GIPM Y (OMNI mean)'] > y_lim_sw[0]) & (df_high['GIPM Y (OMNI mean)'] < y_lim_sw[1])]
    sw_low = df_low.loc[(df_low['GIPM X (OMNI mean)'] > x_lim_sw[0]) & (df_low['GIPM X (OMNI mean)'] < x_lim_sw[1]) & (df_low['GIPM Y (OMNI mean)'] > y_lim_sw[0]) & (df_low['GIPM Y (OMNI mean)'] < y_lim_sw[1])]
    
    para_msh_high = df_high.loc[(df_high['GIPM X (OMNI mean)'] > x_lim_para_msh[0]) & (df_high['GIPM X (OMNI mean)'] < x_lim_para_msh[1]) & (df_high['GIPM Y (OMNI mean)'] > y_lim_para_msh[0]) & (df_high['GIPM Y (OMNI mean)'] < y_lim_para_msh[1])]
    para_msh_low = df_low.loc[(df_low['GIPM X (OMNI mean)'] > x_lim_para_msh[0]) & (df_low['GIPM X (OMNI mean)'] < x_lim_para_msh[1]) & (df_low['GIPM Y (OMNI mean)'] > y_lim_para_msh[0]) & (df_low['GIPM Y (OMNI mean)'] < y_lim_para_msh[1])]
    
    perp_msh_high = df_high.loc[(df_high['GIPM X (OMNI mean)'] > x_lim_perp_msh[0]) & (df_high['GIPM X (OMNI mean)'] < x_lim_perp_msh[1]) & (df_high['GIPM Y (OMNI mean)'] > y_lim_perp_msh[0]) & (df_high['GIPM Y (OMNI mean)'] < y_lim_perp_msh[1])]
    perp_msh_low = df_low.loc[(df_low['GIPM X (OMNI mean)'] > x_lim_perp_msh[0]) & (df_low['GIPM X (OMNI mean)'] < x_lim_perp_msh[1]) & (df_low['GIPM Y (OMNI mean)'] > y_lim_perp_msh[0]) & (df_low['GIPM Y (OMNI mean)'] < y_lim_perp_msh[1])]

    sample_dict = {"Foreshock high":fs_high, "Foreshock low": fs_low, "SW high":sw_high, "SW low": sw_low,"Qpara MSH high": para_msh_high, "Qpara MSH low": para_msh_low,"Qperp MSH high": perp_msh_high, "Qperp MSH low": perp_msh_low}
    
    return(sample_dict)


def power_comparison(sample_dict, quantity):

    #upstream
    
    fs_high = sample_dict["Foreshock high"]
    fs_low = sample_dict["Foreshock low"]
    sw_high = sample_dict["SW high"]
    sw_low = sample_dict["SW low"]

    #m'sheath
    qpara_high = sample_dict["Qpara MSH high"]
    qpara_low = sample_dict["Qpara MSH low"]
    qperp_high = sample_dict["Qperp MSH high"]
    qperp_low = sample_dict["Qperp MSH low"]

    print("\033[1m" + 'Std Dev Normalised Transverse Power by Region, alpha=30-45deg' + "\033[0;0m") 
    print(f'Foreshock, High {quantity}:', fs_high['ULF Band Normalised Transverse Power'].std())
    print(f'Foreshock, Low {quantity}:', fs_low['ULF Band Normalised Transverse Power'].std())
    print(f'Q-para MSH, High {quantity}:', qpara_high['ULF Band Normalised Transverse Power'].std())
    print(f'Q-para MSH, Low {quantity}:', qpara_low['ULF Band Normalised Transverse Power'].std())
    print(' ')
    print(f'Solar Wind, High {quantity}:', sw_high['ULF Band Normalised Transverse Power'].std())
    print(f'Solar Wind, Low {quantity}:', sw_low['ULF Band Normalised Transverse Power'].std())
    print(f'Q-perp MSH, High {quantity}:', qperp_high['ULF Band Normalised Transverse Power'].std())
    print(f'Q-perp MSH, Low {quantity}:', qperp_low['ULF Band Normalised Transverse Power'].std())
    
    print(' ')
    print("\033[1m" + 'Std Dev Normalised Compressive Power by Region' + "\033[0;0m") 
    print(f'Foreshock, High {quantity}:', fs_high['ULF Band Normalised Compressive Power'].std())
    print(f'Foreshock, Low {quantity}:', fs_low['ULF Band Normalised Compressive Power'].std())
    print(f'Q-para MSH, High {quantity}:', qpara_high['ULF Band Normalised Compressive Power'].std())
    print(f'Q-para MSH, Low {quantity}:', qpara_low['ULF Band Normalised Compressive Power'].std())
    print(' ')
    print(f'Solar Wind, High {quantity}:', sw_high['ULF Band Normalised Compressive Power'].std())
    print(f'Solar Wind, Low {quantity}:', sw_low['ULF Band Normalised Compressive Power'].std())
    print(f'Q-perp MSH, High {quantity}:', qperp_high['ULF Band Normalised Compressive Power'].std())
    print(f'Q-perp MSH, Low {quantity}:', qperp_low['ULF Band Normalised Compressive Power'].std())
    
    print("\033[1m" + 'Mean \u00B1 Standard Error, Transverse Power, alpha=30-45deg' + "\033[0;0m") 
    print(f'Foreshock, High {quantity}:', fs_high['ULF Band Normalised Transverse Power'].mean(), u"\u00B1", fs_high['ULF Band Normalised Transverse Power'].sem())
    print(f'Foreshock, Low {quantity}:', fs_low['ULF Band Normalised Transverse Power'].mean(), u"\u00B1", fs_low['ULF Band Normalised Transverse Power'].sem())
    print(f'Q-para MSH, High {quantity}:', qpara_high['ULF Band Normalised Transverse Power'].mean(), u"\u00B1", qpara_high['ULF Band Normalised Transverse Power'].sem())
    print(f'Q-para MSH, Low {quantity}:', qpara_low['ULF Band Normalised Transverse Power'].mean(), u"\u00B1", qpara_low['ULF Band Normalised Transverse Power'].sem())
    print(' ')
    print(f'Solar Wind, High {quantity}:', sw_high['ULF Band Normalised Transverse Power'].mean(), u"\u00B1", sw_high['ULF Band Normalised Transverse Power'].sem())
    print(f'Solar Wind, Low {quantity}:', sw_low['ULF Band Normalised Transverse Power'].mean(), u"\u00B1", sw_low['ULF Band Normalised Transverse Power'].sem())
    print(f'Q-perp MSH, High {quantity}:', qperp_high['ULF Band Normalised Transverse Power'].mean(), u"\u00B1", qperp_high['ULF Band Normalised Transverse Power'].sem())
    print(f'Q-perp MSH, Low {quantity}:', qperp_low['ULF Band Normalised Transverse Power'].mean(), u"\u00B1", qperp_low['ULF Band Normalised Transverse Power'].sem())
    
    print(' ')
    print("\033[1m" + 'Mean \u00B1 Standard Error, Compressive Power, alpha=30-45deg' + "\033[0;0m") 
    print(f'Foreshock, High {quantity}:', fs_high['ULF Band Normalised Compressive Power'].mean(), u"\u00B1", fs_high['ULF Band Normalised Compressive Power'].sem())
    print(f'Foreshock, Low {quantity}:', fs_low['ULF Band Normalised Compressive Power'].mean(), u"\u00B1", fs_low['ULF Band Normalised Compressive Power'].sem())
    print(f'Q-para MSH, High {quantity}:', qpara_high['ULF Band Normalised Compressive Power'].mean(), u"\u00B1", qpara_high['ULF Band Normalised Compressive Power'].sem())
    print(f'Q-para MSH, Low {quantity}:', qpara_low['ULF Band Normalised Compressive Power'].mean(), u"\u00B1", qpara_low['ULF Band Normalised Compressive Power'].sem())
    print(' ')
    print(f'Solar Wind, High {quantity}:', sw_high['ULF Band Normalised Compressive Power'].mean(), u"\u00B1", sw_high['ULF Band Normalised Compressive Power'].sem())
    print(f'Solar Wind, Low {quantity}:', sw_low['ULF Band Normalised Compressive Power'].mean(), u"\u00B1", sw_low['ULF Band Normalised Compressive Power'].sem())
    print(f'Q-perp MSH, High {quantity}:', qperp_high['ULF Band Normalised Compressive Power'].mean(), u"\u00B1", qperp_high['ULF Band Normalised Compressive Power'].sem())
    print(f'Q-perp MSH, Low {quantity}:', qperp_low['ULF Band Normalised Compressive Power'].mean(), u"\u00B1", qperp_low['ULF Band Normalised Compressive Power'].sem())

    quantity_dict = {'MA': 'M_A (mean)', 'B':'IMF B (mean)', 'V': 'SW V (mean)', 'Np': 'SW Np (mean)'}

    quant_ex = quantity_dict[quantity]
    
    print("\033[1m" + 'Foreshock Changes:' + "\033[0;0m") 
    print(f'Mean {quantity}, high:', fs_high[quant_ex].mean(),f'Mean Normalised Transverse Power, high {quantity}:', fs_high['ULF Band Normalised Transverse Power'].mean() ,f'Mean Compressive, high {quantity}:', fs_high['ULF Band Normalised Compressive Power'].mean())
    print(f'Mean {quantity}, low:', fs_low[quant_ex].mean(),f'Mean Normalised Transverse Power, low {quantity}:', fs_low['ULF Band Normalised Transverse Power'].mean() ,f'Mean Compressive, low {quantity}:', fs_low['ULF Band Normalised Compressive Power'].mean())
    print(f'Ratio, {quantity} Change:', fs_high[quant_ex].mean()/fs_low[quant_ex].mean())
    print('Ratio, Normalised Transverse Change:', fs_high['ULF Band Normalised Transverse Power'].mean()/fs_low['ULF Band Normalised Transverse Power'].mean())
    print('Ratio, Normalised Compressive Change:', fs_high['ULF Band Normalised Compressive Power'].mean()/fs_low['ULF Band Normalised Compressive Power'].mean())
    
    
    print("\033[1m" + 'Foreshock vs Solar Wind Compressibility' + "\033[0;0m") 
    print(f'Mean FS Compressibility, Low {quantity}:', fs_low['Compressibility'].mean(),f'Mean SW Compressibility, Low {quantity}:', sw_low['Compressibility'].mean() ,'Ratio:', (fs_low['Compressibility'].mean())/(sw_low['Compressibility'].mean()))

    return()