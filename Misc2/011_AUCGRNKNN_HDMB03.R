## This script obtains NMF representation od PDX tumor cells using NMF model
## trained using patient tumor data in TF-GRN AUC value sapce

library(reticulate)

Sys.setenv(RETICULATE_PYTHON = "/path/to/python")
use_python("/path/to/python")

#function for  return nearest neighbour less than radius and max neighbors
suppressPackageStartupMessages({
  library(Matrix)
  library(uwot)
  library(scater)
  library(AUCell)
})
DIR_GRN <- "GRNana/"
setwd(DIR_GRN)
#
## intial processing ----
load(paste0(DIR_GRN,"pdG34MB_G34GRN_AUCv2.RData") )
pd_ref <- plot.data.com
ms_ref <- t(assay(ms))

## just focusing on first pair of HDMB03 PDX
load(paste0(DIR_GRN,"pdHDMB03_G34GRN_AUCv2.RData") )
ms_tar <- t(assay(ms))
ms_tar <- ms_tar[,colnames(ms_ref)]


###python -----
repl_python()
import sklearn.decomposition as sk
import numpy as np
import random

model = sk.NMF(n_components=25,init="nndsvd", max_iter=10000, random_state=0, tol=1e-5) 
random.seed(456)
model.fit(r.ms_ref)
wn = model.transform(r.ms_tar)

exit
rd <- py$wn
set.seed(134)
del1 <- umap(rd, n_components = 2, min_dist = 0.1, metric = "cosine", n_neighbors = 25)
colnames(del1) <-  c("nfUMAP1","nfUMAP2")
set.seed(134)
del2 <- umap(rd, min_dist = 0.1, n_neighbors = 25)
colnames(del2) <-  c("nfUMAP1b","nfUMAP2b")
plot.data <- cbind(plot.data,del1,del2)
##knn ##------

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

random.seed(455)
a = sk.kneighbors_graph(wn,n_neighbors=26, metric="cosine", n_jobs=3, include_self=True) 
nxg = nx.from_scipy_sparse_matrix(a) 
g_rnn = ig.Graph(len(nxg), list(zip(*list(zip(*nx.to_edgelist(nxg)))[:2])), directed=False)
del(nxg,a)

ela = g_rnn.get_edgelist()
wj = g_rnn.similarity_jaccard(pairs=ela) #for edge list pairs
part = la.find_partition(g_rnn, la.ModularityVertexPartition, weights=wj, n_iterations= -1, max_comm_size=1000, seed=12)

k6 = part.modularity
#clusters
k7 = np.array(part.membership)

exit
k10 <- py$k7

print(py$k6)
#assigning clustering-
plot.data$KNN_cl <- paste0("msKNN_",k10)

save(plot.data, file = paste0(DIR_GRN,"pdHDMB03_G34MBaucKNN.RData"))

q()
