import os
import sys
import json

# from numba import jit

import numpy as np
import pandas as pd

# import matplotlib.pyplot as plt
# import matplotlib.colors as colors 
# from matplotlib import cm

# from tqdm import tqdm

#===================================================================================== 
def read_gauss_params_data(filepath):
    """

    assumes format: 
    """
    gauss_paramsets = []
    with open(filepath) as f:
        for line in f.readlines():
            if len(line) > 2:
                gauss_paramsets.append(json.loads(line[:-2].replace("'",'"')))
    
    mus, ns, Ms = [], [], []
    for row in gauss_paramsets:
        mus.append(row["mu"])
        ns.append(row["n"])
        Ms.append(np.linalg.inv(row["sigma"]))
    
    return mus, ns, Ms

def decode_filename(filename):
    """
    assumes format
    """
    aa_info, num_of_gaussians, _ = filename.split(".")
    aa_code = aa_info[0]
    ifpreproline = (aa_info[1].lower() == "p")
    return aa_code, ifpreproline, int(num_of_gaussians)

#=====================================================================================
def create_grid_of_function_values(input_function, grid_res=256, fun_arg_range=[-np.pi, np.pi]):
    """
    """
    ls = np.linspace(fun_arg_range[0], fun_arg_range[1], grid_res)
    xx,yy = np.meshgrid(ls, ls)
    samp_points = np.transpose(np.array([xx.flatten(), yy.flatten()]))
    samp_probs = input_function(samp_points)
    arr = samp_probs.reshape((grid_res,-1))
    return arr

def get_normalization_factor(input_function, fun_arg_range=[-np.pi, np.pi], test_grid_res=256):
    """
    """
    test_grid_res = 256
    dx = 2 * np.pi / test_grid_res   
    arr = create_grid_of_function_values(input_function, grid_res=256)
    vol = np.sum(arr) * dx ** 2
    return vol

# @jit
def get_probs_nb(points, nns, mus, Ms, n, l):
    """
    """
    probs = np.zeros(n)
    for ip in range(n):
        for ig in range(l):
            d = points[ip] - mus[ig]
            probs[ip] += nns[ig] * np.exp(-(d @ Ms[ig] @ d))
    return probs



