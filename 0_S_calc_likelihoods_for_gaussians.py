import os
import sys
import numpy as np
import pandas as pd

import concurrent
from concurrent.futures import ProcessPoolExecutor

from tqdm import tqdm

# from my_utilities.constants import AMINOACIDS
from my_utilities.gaussians import *
from my_utilities.data import get_phipsi_dict_from_dataframe

# READ COMMAND LINE ARGUMENTS
#======================================================================================================================#

if len(sys.argv) != 4:
    sys.exit(f"Usage: {sys.argv[0]} [directory path] [aa ramachandran data .pkl] [output path (.pkl)]")

gauss_params_dir = sys.argv[1]
aa_phipsi_data = get_phipsi_dict_from_dataframe(pd.read_pickle(sys.argv[2]))
output_path = sys.argv[3]

# FOR ALL THE FILES WE CALCULATE THE LOG-LIKELIHOODS
#======================================================================================================================#

def calculate_likelihood_for_gmm(fname):
    """
    """
    print(f"Starting calculations for {fname}")

    aa_code, ifpreproline, num_of_gaussians = decode_filename(fname)
    try:
        mus, ns, Ms = read_gauss_params_data(f"{gauss_params_dir}/{fname}")
    except Exception as e:
        print(f"WARNING: Encountered issue, while trying to open file {fname}")
        print(e, file=sys.stderr, flush=True)
        return None

    if len(ns) != num_of_gaussians:
        print(f"WARNING: Wrong input file lengh; file: {fname}, length: {len(ns)}")
        return None
    
    else:
        prob_density_not_norm = lambda points : get_probs_nb(points, ns, mus, Ms, len(points), len(ns))
        
        # NORMALIZATION
        normalization_factor = get_normalization_factor(prob_density_not_norm, fun_arg_range=[-np.pi, np.pi], test_grid_res=256)
        if normalization_factor == 0:
            sys.exit(f"Error: zero volume: {fname}")

        nns = [n / normalization_factor for n in ns] 
        prob_density = lambda points : get_probs_nb(points, nns, mus, Ms, len(points), len(ns))
        
        # GET ACTUALLY OBSERVED POINTS (PHI,PSI)
        observed_points = aa_phipsi_data[aa_code, ifpreproline]

        # CALCULATE LIKELIHOOD
        log_likelihood = np.sum(np.log(prob_density(observed_points)))

        # SAVE DATA
        num_of_observations = observed_points.shape[0] # number of observed points
        norm_log_likelihood = log_likelihood / num_of_observations
        n_params = 6 * num_of_gaussians # amplitude, center (2D), 2x2 sigma matrix, but symmetrical, so 3 params

        print(f"Finished calculations for {fname}")

        return {"aminoacid": aa_code,
                "preproline": ifpreproline,
                "representation": "gaussians",
                "num_of_gaussians": num_of_gaussians,
                "log_likelihood": log_likelihood,
                "norm_log_likelihood": norm_log_likelihood,
                "parameters": n_params,
                "data points": num_of_observations}


filenames = os.listdir(gauss_params_dir)

log_likelihood_data_list = []

with ProcessPoolExecutor() as executor:
    futures = [executor.submit(calculate_likelihood_for_gmm, filename) for filename in filenames]
 
    for future in concurrent.futures.as_completed(futures):
        log_likelihood_data_list.append(future.result())

log_likelihood_data_list_corr = [llinfo for llinfo in log_likelihood_data_list if llinfo is not None]

# WE DUMP THE RESULTS INTO DATAFRAME AND INTO A PICKLE
log_likelihood_data = pd.DataFrame(log_likelihood_data_list_corr)
log_likelihood_data.to_pickle(output_path)