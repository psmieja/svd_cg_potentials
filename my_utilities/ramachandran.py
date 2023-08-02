import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm

def get_ramachandran(data, aminoacid_code,
                          ifpreproline=False,
                          n_bins=32,
                          excluded_volume=None, # excluded volume mask
                          norm=True,
                          pseudocount=0.001):
    """
    """
    dx = 2 * np.pi / n_bins
    bins = [-np.pi + i*dx for i in range(n_bins+1)]
    
    our_data = data[(data["aa_1l"]==aminoacid_code) & (data["preproline"] == ifpreproline)][["phi","psi"]]
    ramachandran = np.histogram2d(our_data.values[:,0], our_data.values[:,1], bins=bins)
    hist = ramachandran[0]

    if excluded_volume is not None:
        hist = hist * excluded_volume
        
    hist += pseudocount
    
    # histogram normalization ==> volume should be equal to 1
    if norm:
        hist = hist / (np.sum(hist) * dx * dx)
    
    return hist, n_bins, dx, our_data

def plot_ramachandran(arr, save_to=None):
    """
    takes 2D array as input and plots it...
    if save_to is None (default) plot is just shown
    if save_to is a str, plot is saved to its value
    """
    plt.imshow(np.transpose(arr),
               origin="lower",
               norm=colors.LogNorm(vmin=arr.min(), vmax=arr.max()),
               cmap=cm.coolwarm,
               interpolation="none")
    n = arr.shape[0]
    plt.xticks(ticks=[-0.5,n/2-0.5,n-0.5], labels=["$-\pi$", "0", "$\pi$"])
    plt.yticks(ticks=[-0.5,n/2-0.5,n-0.5], labels=["$-\pi$", "0", "$\pi$"])
#     plt.colorbar()
    if save_to is not None:
        plt.savefig(save_to)
    else:
        plt.show()
    plt.close()






