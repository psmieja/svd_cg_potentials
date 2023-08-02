import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm


def plot_ramachandran(hd):
    """
    takes 2D array as input and plots it...
    """
    plt.imshow(np.transpose(hd),
               origin="lower",
               norm=colors.LogNorm(vmin=hd.min(), vmax=hd.max()),
               cmap=cm.coolwarm,
               interpolation="none")
    n = hd.shape[0]
    plt.xticks(ticks=[-0.5,n/2-0.5,n-0.5], labels=["$-\pi$", "0", "$\pi$"])
    plt.yticks(ticks=[-0.5,n/2-0.5,n-0.5], labels=["$-\pi$", "0", "$\pi$"])
#     plt.colorbar()
    plt.show()
    plt.close()