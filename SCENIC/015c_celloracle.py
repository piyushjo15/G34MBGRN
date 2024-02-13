### CellOracle processing 
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import pickle
import celloracle as co


Sample = "MSM144" ## MSM109 or MSM057
print("Start Analysis for Sample: " + Sample)

#First define directories
projDir = os.path.join('SCENIC/GRNind/')
dataDir = os.path.join('SCENIC/RNA/',Sample+'/')
tmpDir = "TMP/"
os.chdir(projDir)

save_folder = "figures_CO"
#load AnData object
adata = sc.read_h5ad(dataDir+'adata_CO.h5ad')
adata

##base GRN from dictionary
with open(dataDir+'TFGRN.pkl', 'rb') as file:
    TFGRN = pickle.load(file)
    
# Instantiate Oracle object
oracle = co.Oracle()
adata.X = adata.raw.X.copy()

# Instantiate Oracle object.
oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name="Annotation",
                                   embedding_name="X_umap")
                                   
# We invert the dictionary above using a utility function in celloracle.
TGtoTFGRN = co.utility.inverse_dictionary(TFGRN)
oracle.addTFinfo_dictionary(TGtoTFGRN)

# Perform PCA
oracle.perform_PCA()

# Select important PCs
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.show()
print(n_comps)
n_comps = min(n_comps, 50)

n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")

k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")

oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)
                      
#cluster-specific GRN for all clusters.
links = oracle.get_links(cluster_name_for_GRN_unit="Annotation", alpha=1,
                         verbose_level=10)
links.filter_links(p=0.001, weight="coef_abs", threshold_number=10000)
# Save Links object.
links.to_hdf5(dataDir+"links.celloracle.links")

oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=1, use_cluster_specific_TFdict=True)

oracle.to_hdf5(dataDir+"output.celloracle.oracle")



