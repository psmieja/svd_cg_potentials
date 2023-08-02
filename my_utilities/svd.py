import numpy as np
import pandas as pd


def SVD_representation(arr, d):
    """
    
    """
    n = arr.shape[0]
    u, s, vh = np.linalg.svd(arr)
    smat_lim = np.zeros((n,n))
    smat_lim[:d, :d] = np.diag(s[:d])
    return u @ smat_lim @ vh