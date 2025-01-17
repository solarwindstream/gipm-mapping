#WIP
import numpy as np
from math import cos, sin, pi
from new_xyz import new_xyz
from merka05_surface_eq_array_GIPM import merka05_surface_eq_array_GIPM

def merka05_closest_point(XMA, r_gipm):

    #fit bow shock
    A = merka05_surface_eq_array_GIPM(XMA)
    
    Xn = r_gipm[0]
    Yn = r_gipm[1]
    Zn = r_gipm[2]
        
    #Find the closest point, Yn and Zn (in GPE) fixed

    DS = 1
    #NOW BEGIN ITERATIONS. EACH ITERATION BEGINS WITH FINDING F AND ITS DERIVATIVES:
    while  (DS > 1e-5):   
        F = A[0]*Xn**2 + A[1]*Yn**2 + A[2]*Zn**2+ 2*A[3]*Xn*Yn + 2*A[4]*Yn*Zn + 2*A[5]*Xn*Zn + 2*A[6]*Xn+2*A[7]*Yn + 2*A[8]*Zn + A[9]
        FX = 2*A[0]*Xn + 2*A[3]*Yn + 2*A[6]
        FY = 2*A[1]*Yn + 2*A[3]*Xn + 2*A[7]
        FZ = 2*A[2]*Zn

        GRADF2=FX**2+FY**2+FZ**2
        DX=FX*F/GRADF2
        DY=FY*F/GRADF2
        DZ=FZ*F/GRADF2
        DS=np.sqrt(DX**2+DY**2+DZ**2)

        Xn=Xn-DX
        Yn=Yn-DY
        Zn=Zn-DZ

    r_bs_gipm = np.array([Xn, Yn, Zn])
    
    x = r_bs_gipm[0]
    y = r_bs_gipm[1]
    z = r_bs_gipm[2]
    normal_gipm = np.array([2*x*A[0]+2*y*A[3]+2*z*A[5]+2*A[6],2*x*A[3]+2*y*A[1]+2*z*A[4]+2*A[7],2*x*A[5]+2*y*A[4]+2*z*A[2]+2*A[8]])
    normal_gipm = normal_gipm/np.linalg.norm(normal_gipm)

    return r_bs_gipm, normal_gipm