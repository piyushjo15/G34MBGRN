#creating andata object for analysis with CellOracle
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")


suppressPackageStartupMessages({
  library(scater)
  library(scran)
})


#load RNA SCE
Sample <- commandSample(TRUE)
msn <- gsub("MSM","MSN",Sample)
##Dir
DIR_RNA <- "RNA/"
DIR_GRNind <- "SCENIC/GRNind/"

##load RNA DATA
seRNA <- get(load(paste0(DIR_RNA,msn,"scepro.RData")))

#load plot.data (combined RNA) and hvg list
com.vg <- readLines("Files/topHVGG34new10xTUMnoSX_1500.txt") 
##load UMAP from TF-GRN
load(paste0(DIR_GRNind,msn,"_AUC_DM_UMAP.RData"))

#subset for tumor cells of selected sample
cellMeta <- pd[!is.na(pd$ind_cluster),]
table(colnames(seRNA) %in% row.names(cellMeta))
seRNA <- seRNA[com.hvg,row.names(cellMeta)]
#convert RNA SCE to Scanpy AnnData
############# step-by-step to Anndata ###############

#save the matrix
a <- seRNA@assays@data$counts
colnames(a) <- colnames(seRNA)

#save the barcodes
barcodes<-data.frame(colnames(seRNA))
colnames(barcodes)<-'Barcode'

#Save the gene names
genes<-data.frame(rownames(seRNA))
colnames(genes)<-'Gene'

##pca and umap
pca <- reducedDim(seRNA,"PCA")
umap <- as.matrix(cellMeta[,c("aucUMAP1","aucUMAP2")])
dm <- as.matrix(cellMeta[,c("DC1","DC2")])

######### Save the data to python ############
repl_python()

import scanpy as sc
import pandas as pd
from scipy import io
import os

Sample = r.Sample
projDir = os.path.join("SCENIC/RNA/",Sample)

#log count Matrix 
counts = r.a
from scipy import sparse
counts=sparse.csr_matrix(counts)

#load metadata
barcodes=r.barcodes
genes=r.genes

adata=sc.AnnData(counts.T)
adata.raw = adata
adata.obs_names=barcodes['Barcode'].values
adata.var_names=genes['Gene'].values

#import metadata
adata.obs=r.cellMeta
adata.obsm['X_umap']=r.umap
adata.obsm['X_pca']=r.pca
adata.obsm['X_dmap']=r.dm
#save adata object
adata.write(os.path.join(projDir, 'adata_CO.h5ad'), compression='gzip')
exit
print(paste0("Saved adata object for sample:", Sample))
q()
