## This script creates anData for Retina RNA data

#create AnnData Non-Multiome Sample
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")

suppressPackageStartupMessages({
  library(ArchR)
  library(tibble)
  library(hexbin)
  library(SingleCellExperiment)
  library(BSgenome.Hsapiens.UCSC.hg38)
})


#load RNA SCE
Sample <- commandArgs(TRUE)

# Define the folder path and name
projDir = "Retina/SCENIC/"
tmpDir = 'TMP/'
DIR_RNA = "Retina/RNA/"
setwd(projDir)
#load RNA data
seRNA <- get(load(paste0(DIR_RNA,Sample,"scepro.RData")))
##fixing cell Ids
colnames(seRNA) <- paste0(Sample,"_",colnames(seRNA))
#load hvg list
load(paste0(DIR_RNA,"topHVGRetRiboMTSX_PC.RData")) #topHVG
load(paste0(DIR_RNA,"plotdataRetANN.RData"))
#remove ND cells
keep <- plot.data$Annotation=="ND"
plot.data <- plot.data[!keep,]
#convert RNA SCE to Scanpy AnnData
############ add missing information and filter for ATAC cells ###############

plot.data_RNA <- plot.data[plot.data$Batch==Sample,]


######## subset RNA genes to contain only HVGs ############
#filter for HVG
#subst the SCE to contain only the top HVG and the cells also present in the ATAC project
seRNA <- seRNA[rownames(seRNA) %in% topHVG,row.names(plot.data_RNA)]

############ cp the ind_cluster as a new column naming predictedGroup (to have the same name as in ATAC data)

plot.data_RNA$predictedGroup <- plot.data_RNA$Annotationlv2


############# step-by-step to Anndata ###############

#save the matrix
a <- seRNA@assays@data$logcounts
colnames(a) <- colnames(seRNA)

#save the barcodes
barcodes<-data.frame(colnames(seRNA))
colnames(barcodes)<-'Barcode'

#Save the gene names
genes<-data.frame(rownames(seRNA))
colnames(genes)<-'Gene'

#save metadata
cellMeta<-plot.data_RNA

######### Save the data to python ############
repl_python()

import scanpy as sc
import pandas as pd
from scipy import io
import os

Sample = r.Sample
projDir = os.path.join("Retina/SCENIC/",Sample)

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
cellMeta=r.cellMeta
adata.obs=cellMeta
adata.obsm['X_umap']=adata.obs.loc[:,['iUMAP1','iUMAP2']].values.copy() 

#save adata object
adata.write(os.path.join(projDir, 'adata.h5ad'), compression='gzip')

print("Saved adata object for sample: " + Sample)
exit
q()
