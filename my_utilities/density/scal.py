from scipy.interpolate import RegularGridInterpolator
import numpy as np



def get_g_hist_trivial(hist_norm, dx):
    """
    returns probability density function based on simple approximation from grid
    point is some iterable containing phi and psi
    """
    def density_function(point):
        i, j = np.floor((point - np.pi)/dx)
        return hist_norm[int(i)][int(j)]
    
    return density_function



def get_g_hist_linear(hist_norm, dx, n_bins):
    """
    returns probability density function based on simple approximation from grid
    point is some iterable containing phi and psi
    """
    
    # first we use periodic boundary conditions to extend the array
    # so that we can interpolate for the values close to the original array's edge
    
    hist_ext1 = np.vstack([hist_norm[-1,:][None,:],
                       hist_norm,
                       hist_norm[0,:][None,:]])

    hist_ext2 = np.hstack([hist_ext1[:,-1][:,None],
                       hist_ext1,
                       hist_ext1[:,0][:,None]])
    
    #
    
    ext_bin_values = [-np.pi + (i - 0.5)*dx for i in range(n_bins+2)]
    
    
    interp = RegularGridInterpolator((ext_bin_values, ext_bin_values), hist_ext2)
    #...Interpolator is meant to be called on arrays and we want to use it on single points
    interp1 = lambda point: interp(point)[0] 
    return interp1