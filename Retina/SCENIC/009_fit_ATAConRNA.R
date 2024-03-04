## This script calculates projection of ATAC cells on RNA cells using predicted
## AUC scores
library(reticulate)

Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")

#function for  return nearest neighbour less than radius and max neighbors
suppressPackageStartupMessages({
  library(Matrix)
  library(uwot)
  library(batchelor)
  library(BiocParallel)
})
DIR_RNA <- "Retina/RNA/"
DIR_GRN <- "Retina/ATAC/SCENIC/scenicplus/"
DIR_ATAC <- "Retina/ATAC/"

## intial processing ----

load(paste0(DIR_GRN,"pdRetR_RetaucGRN_AUC.RData"))
rm(plot.data)
ms_ref <- scr
rm(auc, scr, plot.data)
#ATAC
load(paste0(DIR_GRN,"imputed_aucGRN_ATACn5.RData"))
rm(pd_ATAC)

load(paste0(DIR_ATAC, "com_RNA_ATAC_pd_scJ.RData"))


###python -----
repl_python()
import sklearn.decomposition as sk
import numpy as np
import random

model = sk.NMF(n_components=25,init="nndsvd", max_iter=10000, random_state=0, tol=1e-5)
random.seed(456)
model.fit(r.ms_ref)
wn = model.transform(r.ms_ref)
wn2 = model.transform(r.ms_atac)

exit
rd_ref <- py$wn
rd_tar <- py$wn2
rd <- rbind(rd_ref,rd_tar)

set.seed(134)
mnn <- reducedMNN(rd, batch=plot.data$Platform, BPPARAM= MulticoreParam(6))

rd_cor <- mnn$corrected

set.seed(134)
del1 <- umap(rd, min_dist = 0.8, metric = "cosine", n_neighbors = 25)
colnames(del1) <-  c("nfUMAP1","nfUMAP2")
plot.data <- cbind(plot.data,del1)

save(rd,rd_cor, plot.data, file =paste0(DIR_GRN,"pdRNA_ATAC_fMNN_AUC.RData"))
q()
