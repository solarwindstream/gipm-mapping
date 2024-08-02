import numpy as np
from math import pi, cos, sin

def merka05_surface_eq_array_GIPM(XMA):

    #Parameters and conventions
    #D solar wind proton number density #/cm3
    #V solar wind speed km/s
    #XMA Alfvenic Mach number

    #Calculating relevant upstream values
    #Vabs = V
    #D = n

    #Merka05 model parameters
    #---------------------------
    #Scaling factor of the coordinate system
    #FAC=((D/7.0)*(Vabs/457.5)**2)**0.166667

    #Fitting parameters for Alfven Mach number dependence
    B11 =0.1089
    B12 =-0.0146
    B31 =0.8063
    B32 =0.0203
    B41 =-0.0302
    B42 =0.0032
    B71 =14.35
    B72 =0.4555
    B73 =111.2
    B81 =0.6111
    B82 =-0.0397
    B101 =-343.6
    B102 =-12.306
    B103 =-3290

    #Fitting parameters of the shock surface
    A1=B11+B12*XMA
    A2=1.0*np.ones(np.shape(XMA))            #As is Peredo
    A3=B31+B32*XMA
    A4=B41+B42*XMA
    A5=0.0*np.ones(np.shape(XMA))            #Assume north-south symmetry -> A5 = A6 = A9 = 0
    A6=0.0*np.ones(np.shape(XMA))
    A7=B71+B72*XMA+B73/((XMA-1)**2)
    A8=B81+B82*XMA
    A9=0.0*np.ones(np.shape(XMA))
    A10=B101+B102*XMA+B103/((XMA-1)**2)

    A = np.array([A1, A2, A3, A4, A5, A6, A7, A8, A9, A10])

    return A

