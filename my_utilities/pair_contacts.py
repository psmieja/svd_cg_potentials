import os
# from numba import jit
from tqdm import tqdm

import numpy as np
import pandas as pd

import sklearn 
from sklearn.metrics import matthews_corrcoef

from scipy.optimize import minimize

import matplotlib.pyplot as plt
import matplotlib.colors as colors 
from matplotlib import cm

from my_utilities.constants import AMINOACIDS, aminoacid_name_from_1l

def load_data(dir_path):
    """
    load all .dat files containing contact data from given dir and return in a DataFrame
    """
    fnames = [f for f in os.listdir(dir_path) if f.split(".")[-1] == "dat"]
    dfs = []
    for fname in tqdm(fnames):
        try:
            dfs.append(
                pd.read_csv(dir_path + fname,
                            delim_whitespace=True,
                            index_col=False,
                            names=["aa1", "aa2", "dij", "dcg", "daa"],
                            skiprows=1)
            )
        except:
            print(f"Problem with file: {fname}")   
    return pd.concat(dfs)

def plot_contact_hists(dfx, N_BINS=100, plot_mcc_curve=False, save_to=None):
    """
    """
    aac_hist, aac_edges = np.histogram(dfx[dfx["aa_contact"]]["dcg"], N_BINS)
    aac_rdf = aac_hist / ((3/4)*(np.pi)*aac_edges[1:]**3)

    aanc_hist, aanc_edges = np.histogram(dfx[~dfx["aa_contact"]]["dcg"], N_BINS)
    aanc_rdf = aanc_hist / ((3/4)*(np.pi)*aanc_edges[1:]**3)

    fig, ax = plt.subplots()
    ax.plot(aac_edges[1:], aac_rdf, label="contact")
    ax.plot(aanc_edges[1:], aanc_rdf, label="no contact")
    ax.tick_params(direction="in")
    ax.set_xlabel("d [Å]")
    ax.set_ylabel("n")
    ax.legend()

    if plot_mcc_curve:
        rs, mccs = [], []
        for r in tqdm(np.linspace(min(aac_edges[0], aanc_edges[0]),
                                  max(aac_edges[-1], aanc_edges[-1]), 100)):
            rs.append(r)
            mccs.append(
                sklearn.metrics.matthews_corrcoef(dfx["aa_contact"], dfx["dcg"] < r)
            )
        ax1 = ax.twinx()
        ax1.plot(rs, mccs, color="red",linestyle="dashed")
        ax1.tick_params(axis='y', labelcolor="red")
        ax1.set_ylabel('MCC', color="red")

    if save_to is None:
        plt.show()
    else:
        plt.savefig(save_to)
    
def get_optimal_d2(dfx):
    """
    function to get the optimal d2 (max length for contact of side-group united atoms)
    such that most effectively predicts all-atom contact
    d2 is basically a parameter of the binary classifier
    and we measure its quality using matthews correlation coefficient (phi)
    """
    realvals = dfx["aa_contact"].to_numpy() # actualy informs whether groups are in contact
    dcgs = dfx["dcg"].to_numpy() # distance between coarse-grained united atoms
    def _fun(r,realvals,dcgs): # we want to maximize the mmc, so we the define the minimized function as its negative
        return -sklearn.metrics.matthews_corrcoef(realvals, dcgs < r)
    minres = minimize(_fun,5,args=(realvals,dcgs), tol=0.005,method="Powell")
    return minres.x[0]




    # DIR_PATH = "./cg_contacts/"
    # DAA_CONTACT_CUTOFF = 5

    # df = load_data(DIR_PATH)
    # df["aa_contact"] = df.apply(lambda row: row["daa"] < DAA_CONTACT_CUTOFF, axis=1)
    # dfg = df.groupby(["aa1","aa2"])
    # dfnrg = {}
    # for g in dfg: #dfgg.keys():
    #     g_key, g_df = g
    #     if (g_key[1],g_key[0]) not in dfnrg.keys():
    #         dfnrg[g_key] = pd.concat([dfg.get_group(g_key), dfg.get_group((g_key[1],g_key[0]))])

    # for
    

