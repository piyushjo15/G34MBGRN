##This script performs in silico gene knock-out in a specific sample
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import pickle
import celloracle as co

Sample = sys.argv[1] #MSM057, MSM109, MSM144

print("Start Analysis for Sample: " + Sample)

#First define directories
projDir = os.path.join('SCENIC/GRNind/')
dataDir = os.path.join('SCENIC/RNA/',Sample+'/')

tmpDir = "TMP/"
os.chdir(projDir)
save_folder = "figures_CO"

oracle = co.load_hdf5(dataDir+"output.celloracle.oracle")
links = co.load_hdf5(dataDir+"links.celloracle.links")

goi = "CRX" # or EOMES
oracle.simulate_shift(perturb_condition={goi: 0.0}, n_propagation=3)

# Get transition probability
oracle.estimate_transition_prob(n_neighbors=200,
                                knn_random=True,
                                sampled_fraction=1)

# Calculate embedding
oracle.calculate_embedding_shift(sigma_corr=0.05)

oracle.to_hdf5(dataDir+"output_"+goi+"KO.celloracle.oracle")
