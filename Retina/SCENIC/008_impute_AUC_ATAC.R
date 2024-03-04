## This script is to impute module score for GRNs using RNA data as a reference.
## I calcualted module score for specific GRNs in the combined Retina RNA data 
## Using scJoint, I obtained the embedding to identify nearest neighbors for each
## ATAC cell in the reference RNA data
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")

#load packages
suppressPackageStartupMessages({
  library(foreach)
  library(scran)
  library(scater)
  library(Matrix)
  library(batchelor)
})
##DIRs
DIR_RNA <- "Retina/RNA/"
DIR_ATAC <- "Retina/ATAC/"
DIR_scJoint <- "Retina/ATAC/scJoint/"
DIR_SCENIC <- "Retina/ATAC/SCENIC/scenicplus/"
##load module score
load(paste0(DIR_SCENIC, "pdRetR_RetaucGRN_AUC.RData"))
rm(plot.data)
ms <- scr
##load combined data
load(paste0(DIR_ATAC, "com_RNA_ATAC_pd_scJ.RData"))
pd_RNA <- plot.data[plot.data$Platform=="RNA",]
pd_ATAC <- plot.data[plot.data$Platform=="ATAC",]
##load embedding from scJoint
RNA_ref <- read.csv(paste0(DIR_scJoint,"output/", "Retina_rna_embeddings.txt"), sep = " ", header = FALSE)
ATAC_tar <- read.csv(paste0(DIR_scJoint,"output/", "Retina_atac_embeddings.txt"), sep = " ", header = FALSE)


################# Imputation ##############################
#using KNN to identify nearest neighbors
## python code ----
repl_python()

import numpy as np
import sklearn.neighbors as skn

neigh = skn.NearestNeighbors(n_neighbors=5) 
neigh.fit(r.RNA_ref)
dis, nd = neigh.kneighbors(r.ATAC_tar, return_distance=True) ##get distance for weighted KNN

exit

##converting into R objects
kneigh <- py$nd ##index of neighbors
kneigh <- kneigh + 1 ##add 1 to index, python vs R differencce
kneigh_dist <- py$dis ##index of neighbors
row.names(kneigh) <- row.names(kneigh_dist) <- row.names(pd_ATAC)
colnames(kneigh) <- colnames(kneigh_dist) <-paste0("N_",seq(dim(kneigh)[2]))

head(kneigh)
head(kneigh_dist)

##converting distance matrix for weight matrix
max_d <- apply(kneigh_dist, 1, max)
kneigh_dist_i <- 1/kneigh_dist #inversing the distance to contribution, farther has lower weight
sum_dis <- rowSums(kneigh_dist_i)
kneigh_w <-  kneigh_dist_i/sum_dis##normalizing  
head(kneigh_w)

##finding psuedocounts from neigbors
ms_atac <- as(matrix(0,dim(ATAC_tar)[1],dim(ms)[2]),"dgCMatrix")
cells_atac <- row.names(ms_atac) <- row.names(pd_ATAC)
colnames(ms_atac) <- colnames(ms)
ms_atac[1:4,1:4]

for(x in cells_atac){
  n_x <- kneigh[x,]
  #weighted
  n_w <- kneigh_w[x,]
  mdt_r <-t(t(ms[n_x,])*n_w)##obtain data from neighbors, multiply each sample by weights and trasnpose
  ms_atac[x,] <- apply(mdt_r,2,sum) ##sum all the expression
  rm(n_x,n_w, mdt_r)

}
rm(x)
save(ms_atac, pd_ATAC,file = paste0(DIR_SCENIC,"imputed_aucGRN_ATACn5.RData" ))
q()