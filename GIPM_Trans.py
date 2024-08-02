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

    #print(X_GIPM, Y_GIPM, Z_GIPM)

    #######need transformation matrix!
    GIPM_trans = np.array([X_GIPM[0], Y_GIPM[0], Z_GIPM[0], X_GIPM[1], Y_GIPM[1], Z_GIPM[1], X_GIPM[2], Y_GIPM[2], Z_GIPM[2]]).reshape((3,3))
    
    return(GIPM_trans)

