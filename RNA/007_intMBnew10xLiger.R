## This script uses LIGER to integrate tumor RNA data
library(reticulate)

Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")
#for integrating CB and tumor data on UMAP for lineage
suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(batchelor)
  library(rlist)
  library(Matrix)
  library(uwot)
  library(BiocParallel)
  library(rliger)
})

# ##load data and process ----
Tum <- list()
Tum2 <- list()
load("plotdataG34_SVM.RData")

sams <- unique(plot.data$Batch)
com.hvg <- readLines("Files/topHVGG34new10xnoRPMTSX.txt")
com.sce_cnnc <- c()
for(x in sams){
  sce <- get(load(paste0(x,"scepro.RData")))
  rm(list=ls(pattern="^MSN"))
  del <- logcounts(sce)
  del <- cosineNorm(del[com.hvg,])
  del2 <- counts(sce)[com.hvg,]
  Tum[[x]] <- round(del2)
  Tum2[[x]] <- as.matrix(t(del))
  rm(del,del2,sce)
}
rm(x)

## Liger
com <- createLiger(Tum)

com@scale.data <- Tum2

com <- optimizeALS(com, k = 50, max.iters=100000, rand.seed=123)

sams <- unique(plot.data$Batch)

del <- com@H
rd_als <- c()
for (x in sams) {
  rd_als <- rbind(rd_als,del[[x]])
}

##UMAP of uncorrected data
set.seed(198)
del_umap <- uwot::umap(rd_als,  n_neighbors = 25, min_dist = 0.3)
plot.data$lUMAP1 <- del_umap[,1]
plot.data$lUMAP2 <- del_umap[,2]

## Batch correcting LIGER factors
set.seed(777)
mnn.nmf <- reducedMNN(rd_als, batch=plot.data$Batch, BPPARAM = MulticoreParam(6))
rd_als_cor <- mnn.nmf$corrected

##UMAP of corrected integrated data
set.seed(134)
del_umap <- uwot::umap(rd_als_cor,  n_neighbors = 25, min_dist = 0.3)
plot.data$lcUMAP1 <- del_umap[,1]
plot.data$lcUMAP2 <- del_umap[,2]


##KNN and leiden clustering -------

repl_python()
import numpy as np
import igraph as ig
import leidenalg as la
import networkx as nx
import scipy.sparse as ss
import sklearn.neighbors as sk

##using sklearn neighbors graph
a = sk.kneighbors_graph(r.rd_als_cor,n_neighbors=11, metric="cosine", n_jobs=3, include_self=True) #+1 n_neighbors

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
k5 = la.find_partition(g_rnn, la.ModularityVertexPartition, weights=wj, n_iterations= -1, max_comm_size=4000, seed=12)

k6 = k5.modularity
#clusters
k7 = np.array(k5.membership)
exit
k10 <- py$k7
print(py$k6)
#assigning clustering-
plot.data$KNN_cl <- paste0("KNN_",k10)

#### identify cells that are non-neuronal nromal
plot.data$Nor_SVMpri <- "No"
keep <- plot.data$SVM_prilabKD %in%  c("astrocyte",'mural/endoth',"immune","meningeal",
                                       "oligo_progenitor","erythroid","oligodendrocyte",
                                       "glioblast")
table(keep)
plot.data[keep1,]$Nor_SVMpri <- "Yes"

keep1 <- plot.data$Nor_SVMpri=="Yes"
table(keep1)
keep2 <- plot.data$Nor_IC=="Yes" ## Normal clusters identified in each sample
table(keep2)
keep3 <- plot.data$Nor_KNNCl=="Yes" ## KNNCl cluster exhibiting non-neural cells
table(keep3)

plot.data$Nor <- "No"
keep <- keep1 | keep2 | keep3
table(keep)
plot.data[keep,]$Nor <- "Yes"


save(plot.data, file = "plotdata_G34MB_SVM_Liger.RData" )
q()
