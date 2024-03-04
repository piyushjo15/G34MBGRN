## Performing NMF factorization followed by RNN and clustering
library(reticulate)

Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")

suppressPackageStartupMessages({
  library(batchelor)
  library(scran)
  library(scater)
  library(Matrix)
})

#
## load multibatchnormed cosnormed  merged----
load("cosnorRet.RData")
load("plotdataRet.RData")


dim(com.sce_cnnc)

mdt <- t(com.sce_cnnc[1:1500,]) ##already did topHVG[1:4000] while cosinenorm
genes <- row.names(com.sce_cnnc)
rm(com.sce_cnnc)
repl_python()
import sklearn.decomposition as sk
import numpy
import random
model = sk.NMF(n_components=25, init='nndsvda', max_iter=10000, random_state=0, tol=1e-5)
random.seed(456)
w = model.fit_transform(r.mdt)
h = model.components_
exit

rd_nmf= py$w
H = py$h

## run UMAP
set.seed(787)
del <- uwot::umap(rd_nmf, n_neighbors=50, min_dist=0.8, metric = "cosine")

plot.data$nfUMAP1 <- del[,1]
plot.data$nfUMAP2 <- del[,2]

## python code ----
repl_python()
import numpy as np
import igraph as ig
import leidenalg as la
import networkx as nx
import random
import sklearn.neighbors as skn

##using sklearn neighbors graph
a = skn.kneighbors_graph(r.rd_nmf,n_neighbors=21, metric="cosine", n_jobs=3, include_self=True) #+1 n_neighbors

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
k5 = la.find_partition(g_rnn, la.ModularityVertexPartition, weights=wj, n_iterations= -1, max_comm_size=2000, seed=12)

k6 = k5.modularity
#clusters
k7 = np.array(k5.membership)
exit
k10 <- py$k7
print(py$k6)
#assigning clustering-
plot.data$KNN_cl <- paste0("KNN_",k10)

save(plot.data, file = "plotdataRetANN.RData" ) ## the obtained clusters were annotated based on markers
q()
