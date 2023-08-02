import numpy as np
from scipy.interpolate import RegularGridInterpolator
import numpy as np


def get_trivial_vec_density(hist_norm):
    
    n = hist_norm.shape[0]
    dx = 2 * np.pi / n

    def density(points):
        idxs   = np.floor((points + np.pi)/dx).astype(int)
        return hist_norm[idxs[:,0],idxs[:,1]]
    
    return density



def get_linear_interpolator(hist_norm):

    n_bins = hist_norm.shape[0]
    dx = 2 * np.pi / n_bins

    hist_ext1 = np.vstack([hist_norm[-1,:][None,:],
                       hist_norm,
                       hist_norm[0,:][None,:]])

    hist_ext2 = np.hstack([hist_ext1[:,-1][:,None],
                       hist_ext1,
                       hist_ext1[:,0][:,None]])

    ext_bin_values = [-np.pi + (i - 0.5)*dx for i in range(n_bins+2)]

    return RegularGridInterpolator((ext_bin_values, ext_bin_values), hist_ext2)