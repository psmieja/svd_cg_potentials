import os
import sys

import concurrent
from concurrent.futures import ProcessPoolExecutor

import skimage
import numpy as np
import pandas as pd

# import matplotlib.pyplot as plt
# import matplotlib.colors as colors
# from matplotlib import cm

from tqdm import tqdm

from my_utilities.data import get_data
from my_utilities.ramachandran import get_ramachandran, plot_ramachandran
from my_utilities.svd import SVD_representation
from my_utilities.density.vec import get_trivial_vec_density, get_linear_interpolator
from my_utilities.metrics import log_likelihood
from my_utilities.constants import AMINOACIDS

from my_utilities.data import get_phipsi_dict_from_dataframe

# CONST PARAMS
#=============
MATRIX_SIZES = [16,32,64,128]

# READ COMMAND LINE ARGS
#=======================================================================================================================
if len(sys.argv) != 4:
    sys.exit(f"Usage: {sys.argv[0]} [phiphsi data .pkl path] [excl vol .pkl path] [output path (.pkl)]")

phipsi_data_path = sys.argv[1]
excluded_volume_mask_path = sys.argv[2]
output_path = sys.argv[3]

aa_phipsi_data = pd.read_pickle(phipsi_data_path)

# EXCLUDED VOLUME
#=======================================================================================================================
print("LOG: Finding the excluded volume matrices")
excluded_volume_mask_high_res = ~pd.read_csv(excluded_volume_mask_path).to_numpy()
excluded_volume_mask_high_res_size =  excluded_volume_mask_high_res.shape[0]
excluded_volume_masks = {}
for matrix_size in MATRIX_SIZES:
    if matrix_size > excluded_volume_mask_high_res_size:
        sys.exit("Error: larger excluded volume mask needed")
    downsampling_factor = int(excluded_volume_mask_high_res_size // matrix_size)
    ta = skimage.measure.block_reduce(excluded_volume_mask_high_res, downsampling_factor)
    excluded_volume_masks[matrix_size] = ta > ta.mean()

# FIND & ELIMINATE THE OUTLIERS (aka points inside the excluded volume)
#=======================================================================================================================
print("LOG: Eliminating the outliers")
n = excluded_volume_mask_high_res_size
dx = 2 * np.pi / n
aa_phipsi_data["outlier"] = aa_phipsi_data.apply(lambda row:
                             ~ excluded_volume_mask_high_res[
                                 tuple(((np.array([row["phi"],row["psi"]]) + np.pi) // dx).astype(int))
                             ],
                             axis=1
                            )
aa_phipsi_data = aa_phipsi_data[~aa_phipsi_data["outlier"]]

# PREPROCESSING
#=======================================================================================================================
print("LOG: Data preprocessing")
aa_phipsi_data = aa_phipsi_data[["phi","psi","aa_1l","preproline"]]
aa_phipsi_data = aa_phipsi_data[aa_phipsi_data["aa_1l"].isin(AMINOACIDS)]

# RAW RAMACHANDRAN PLOT DATA
#=======================================================================================================================
print("LOG: Creating ramachandran matrices")
ramachandran_data = {}
for a in tqdm(AMINOACIDS):
    for n_bins in MATRIX_SIZES:
        for ifpreproline in [True, False]:
            ifpreproline_descr = "PP" if ifpreproline else "NP"
            hist, n_bins, dx, our_data  =  get_ramachandran(aa_phipsi_data, a, 
                                                            ifpreproline=ifpreproline,
                                                            n_bins=n_bins,
                                                            pseudocount=0,
                                                            norm=False)
            ramachandran_data[(a, ifpreproline, n_bins)] = {"hist": hist, "raw_data": our_data}

# FIND SVD REPRESENTATIONS
#=======================================================================================================================
print("LOG: Calculating the SVD representations")
svd_representations = {}
for aa in tqdm(AMINOACIDS):
    for ifpp in [True, False]: #ifpreproline
        for n in MATRIX_SIZES:
            dx = 2 * np.pi / n
            rd = ramachandran_data[(aa,ifpp,n)]
            hist = rd["hist"]
            for d in range(int(np.floor(n//2))):
                hist_svd = SVD_representation(hist, d) \
                           * excluded_volume_masks[n] \
                # remove negative values, add pseudocount
                hist_svd = np.where(hist_svd > 0, hist_svd, 0.1)
                # normalize
                hist_svd /= np.sum(hist_svd) * (2 * np.pi / n)**2
                svd_representations[aa,ifpp,n,d] = hist_svd

# CALCULATE LIKELIHOODS
#=======================================================================================================================
print("LOG: Calculating likelihoods")

def calculate_likelihood(aa,pp,n,interp_method,d):
    
    rd = ramachandran_data[(aa, pp, n)]
    raw_hist, raw_data = rd["hist"], rd["raw_data"]
    mydata = raw_data.values
    dx = 2 * np.pi / n
    
    if d is None: # non-svd
        representation = "raw matrix"
        n_params = n * n
        myhist = raw_hist * excluded_volume_masks[n]
        myhist = np.where(myhist > 0, myhist, 0.1) # remove negative values, add pseudocount
        myhist /= np.sum(myhist) * dx**2 # normalize

    else: # svd
        representation = "SVD"
        n_params = 2 * n * d + n
        dx = 2 * np.pi / n
        myhist = svd_representations[aa, pp, n, d]
        
    if interp_method == "linear":
        densities = get_linear_interpolator(myhist)(mydata)
    elif interp_method == "trivial":
        densities = get_trivial_vec_density(myhist)(mydata),

    log_likelihood = np.sum(np.log(densities))
    data_points = mydata.shape[0]
    norm_log_likelihood = log_likelihood / data_points

    return {"aminoacid": aa,
              "preproline": pp,
              "representation": representation,
              "interpolation method": interp_method,
              "full matrix size": n,
              "num of SVs": d,
              "log_likelihood": log_likelihood,
              "norm_log_likelihood": norm_log_likelihood,
              "parameters": n_params,
              "data points": data_points}


PARAMS_CUTOFF = 1024 # max number of params per aminoacid
argsets = []
for aa in AMINOACIDS:
    for pp in [True, False]:
        for n in MATRIX_SIZES:
            for interp_method in ["trivial", "linear"]:
                # raw array
                argsets.append([aa,pp,n,interp_method,None])
                # svd representations
                for d in range(int(np.floor(n//2))):                   
                    n_params = 2 * n * d + n
                    if n_params < PARAMS_CUTOFF:
                        argsets.append([aa,pp,n,interp_method,d])

log_likelihood_data_raw = []

with ProcessPoolExecutor() as executor:
    futures = [executor.submit(calculate_likelihood, *args) for args in argsets]
 
    for future in concurrent.futures.as_completed(futures):
        log_likelihood_data_raw.append(future.result())

log_likelihood_data = pd.DataFrame(log_likelihood_data_raw)

# SAVE TO FILE
#=======================================================================================================================
print("LOG: Saving the results")
log_likelihood_data.to_pickle(output_path)