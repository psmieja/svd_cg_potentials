# INSTRUCTION
#
# !mkdir -p ramachandran_plots
# !mkdir -p ramachandran_plots/PP
# !mkdir -p ramachandran_plots/NP
#  ./thisfile dat_files1 ramachandran_plots 32


import sys
import os

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm

from my_utilities.data import get_data
from my_utilities.ramachandran import get_ramachandran, plot_ramachandran


print(sys.argv)

if len(sys.argv) != 4:
    sys.exit(f"Usage: {sys.argv[0]} [datfiles dir] [plots dir] [n bins]")

DATFILES_DIR = sys.argv[1]
PLOTS_DIR = sys.argv[2]
N_BINS = int(sys.argv[3])

data = get_data(DATFILES_DIR)

aminoacids = ['G', 'Q', 'C', 'S', 'V', 'L', 'F', 'P', 'M', 'R', 'N', 'Y', 'E',
              'A', 'D', 'H', 'T', 'I', 'K', 'W']

for ifpreproline in [True, False]:
    ifpreproline_descr = "PP" if ifpreproline else "NP"
    for a in aminoacids:       
        hist, n_bins, dx, our_data  =  get_ramachandran(data, a, 
                          ifpreproline=ifpreproline,
                          n_bins=N_BINS,
                          pseudocount=0)

        plot_ramachandran(hist, save_to=f"{PLOTS_DIR}/{ifpreproline_descr}/{a}.png")


