import numpy as np

def new_xyz(inp, trans_mat):
    """
    Transform the input field to the new frame.

    Parameters
    ----------
    inp : Input data
    trans_mat : Transformation matrix.

    Returns
    -------
    out : Input in the new frame.
    """

    if ((inp.ndim == 1) and (len(inp) == 3)):
        #out = np.matmul(np.matmul(trans_mat.T, inp), trans_mat)
        out = trans_mat.dot(inp)
    else:
        print("Check array dimensions.")

    return out