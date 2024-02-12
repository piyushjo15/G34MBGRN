##using module scores to perform factorization, get integrated UMAp representation
## and then call clusters
library(reticulate)
##python
Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")
#function for  return nearest neighbour less than radius and max neighbors
suppressPackageStartupMessages({
  library(Matrix)
  library(uwot)
  library(scater)
  library(AUCell)
})
DIR_GRN <- "SCENIC/GRNana/"

## loading the AUC score for the combined single-cell tumor data-------
load(paste0(DIR_GRN,"pdG34MB_G34GRN_AUC.RData") )

ms1 <- t(assay(ms)) ##convert into matrix and transpose
ms <- ms1
###python run for NMF factorization-----
repl_python()
import sklearn.decomposition as sk
import numpy as np
import random

model = sk.NMF(n_components=25,init="nndsvd", max_iter=10000, random_state=0, tol=1e-5)
random.seed(456)
model.fit(r.ms)
wn = model.transform(r.ms)

exit
rd <- py$wn
set.seed(134)
del1 <- umap(rd, n_components = 2, min_dist = 0.4, metric = "cosine", n_neighbors = 25)
colnames(del1) <-  c("nfUMAP1","nfUMAP2")
set.seed(134)
plot.data.com <- cbind(plot.data.com,del1)
save(rd, plot.data.com, file=paste0(DIR_GRN,"pdG34MB_G34GRN_AUC_NMF.RData"))

##python run for KNN-leiden clustering------

repl_python()
import numpy as np
import umap
import random
import igraph as ig
import leidenalg as la
import networkx as nx
import scipy.sparse as ss
import sklearn.neighbors as sk

wn=r.rd

##using NMF factors to identify nearest neighbors using sklearn neighbors graph
## This uses data to give a neighborhood graph
random.seed(455)
a = sk.kneighbors_graph(wn,n_neighbors=26, metric="cosine", n_jobs=3, include_self=True) 

## creating a networkx graph
nxg = nx.from_scipy_sparse_matrix(a) 
#converting to igraph
g_rnn = ig.Graph(len(nxg), list(zip(*list(zip(*nx.to_edgelist(nxg)))[:2])), directed=False)
del(nxg,a)

#getting edgelist
ela = g_rnn.get_edgelist()
## gettign weight list
wj = g_rnn.similarity_jaccard(pairs=ela) #for edge list pairs
#getting clustering
part = la.find_partition(g_rnn, la.ModularityVertexPartition, weights=wj, n_iterations= -1, max_comm_size=4000, seed=12)

# #another approach for resolution
# opt = la.Optimiser()
# opt.consider_empty_community = True
# opt.max_comm_size = 1000
# opt.set_rng_seed(123)
# part = la.RBConfigurationVertexPartition(g_rnn,  resolution_parameter = 0.3, weights=wj)
# opt.optimise_partition(part, n_iterations=-1)

k6 = part.modularity
#clusters
k7 = np.array(part.membership)

exit
k10 <- py$k7

print(py$k6)
#assigning clustering-
plot.data.com$KNN_cl <- paste0("msKNN_",k10)

save(plot.data.com, file = paste0(DIR_GRN,"plotdata_G34MBaucKNN.RData"))
q()

