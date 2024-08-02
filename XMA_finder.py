#find average XMA

def XMA_finder(om_averages):
    from math import pi
    import numpy as np
    mu0 = 4*pi*10**(-7) #permittivity vs/(Am)
    m_p = 1.672631*10**(-27) #proton mass kg

    Babs = om_averages.loc[:,'B_mag'].mean()
    Vabs = om_averages.loc[:,'V_gse'].mean()
    Nave = om_averages.loc[:,'Np'].mean()
    rho = Nave*10**6 #in #/m^3

    VA = Babs*10**(-9)/np.sqrt(mu0*rho*m_p)/1000 #km/s
    XMA = Vabs/VA
    
    return(XMA)
