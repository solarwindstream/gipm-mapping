import numpy as np

#GIPM Co-ord Transform

def GIPM_trans(b_gse, v_gse):
    #determine X direction
    ab_correction = np.array([0,-30,0])

    X_GIPM = (-1)*v_gse + ab_correction
    X_norm = np.linalg.norm(X_GIPM)
    X_GIPM = X_GIPM/(X_norm)

    #Y direction 

    #first, B vector

    B_dot_X = np.dot(b_gse,X_GIPM)

    if B_dot_X > 0:
        Y_GIPM = (-1)*b_gse + B_dot_X*X_GIPM
        norm = np.linalg.norm(Y_GIPM)
        Y_GIPM = Y_GIPM/norm
    else:
        Y_GIPM = b_gse - B_dot_X*X_GIPM
        norm = np.linalg.norm(Y_GIPM)
        Y_GIPM = Y_GIPM/norm        
    
    Z_GIPM = np.cross(X_GIPM, Y_GIPM)
    norm = np.linalg.norm(Z_GIPM)
    Z_GIPM = Z_GIPM/norm  
    
    #print(X_GIPM, Y_GIPM, Z_GIPM)
    if not X_GIPM.shape == (3,):
        print('error XGIPM not correct shape')
    if not Y_GIPM.shape == (3,):
        print('error YGIPM not correct shape')  
    if not Z_GIPM.shape == (3,):
        print('error ZGIPM not correct shape')
    
    #check that the magnetic field is entirely in the XY plane:
    
    z_comp = np.dot(b_gse,Z_GIPM)
    z_comp_two_decimals = (round(z_comp, 2))
    
    if not z_comp_two_decimals == 0:
        print('error in GIPM transform, z_comp:', z_comp)

    return(X_GIPM, Y_GIPM, Z_GIPM)

