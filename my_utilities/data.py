import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import matplotlib.colors as colors
from matplotlib import cm

from my_utilities.constants import AMINOACIDS


def get_data(datfiles_dir):
    """
    create pandas array from .dat files in a given directory

    """
    data_list = []
    for filename in os.listdir(datfiles_dir):
        with open(f"{datfiles_dir}{filename}") as file_content:
            for line in file_content.readlines():
                row = line.split()
                # we want to eliminate rows with lacking columns
                if len(row) == 8:
                    data_list.append([str(row[0]),
                                    str(row[1]),
                                    int(row[2]),
                                    str(row[3]),
                                    str(row[4]),
                                    str(row[5]),
                                    float(row[6]),
                                    float(row[7])])
    data =  pd.DataFrame(data_list,
                         columns=["pdb_code",
                                  "chain_code",
                                  "aa_num_in_sequence",
                                  "aa_3l_code",
                                  "aa_triplet",
                                  "2nd_struct",
                                  "phi", "psi"])
    
    # We assign a letter describing the central aminoacid to each row
    data["aa_1l"] = data.apply(lambda row: row["aa_triplet"][1], axis=1)

    # We create a separate column describing whether the aminoacid in that position is preproline
    data["preproline"] = data.apply(lambda row: row["aa_triplet"][2] == "P", axis=1)

    return data

def get_phipsi_dict_from_dataframe(df):
    # from
    res = {}
    for aa in AMINOACIDS:
        for ifpp in [True, False]:
            res[aa,ifpp] = df[(df["aa_1l"]==aa) & (df["preproline"]==ifpp)][["phi","psi"]].to_numpy()

    return res