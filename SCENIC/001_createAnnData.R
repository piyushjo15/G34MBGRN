## This scripts converts the RNA RData into python anData, that will be used as 
## SCENIC+ input
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/ad.dkfz-heidelberg.de/t716r/.conda/envs/scenicplus/bin/python3.8")
use_python("/home/ad.dkfz-heidelberg.de/t716r/.conda/envs/scenicplus/bin/python3.8")

suppressPackageStartupMessages({
  library(ArchR)
  library(tidyr)
  library(dplyr)
  library(SingleCellExperiment)
  library(BSgenome.Hsapiens.UCSC.hg38)
})


#load RNA SCE
args <- commandArgs(TRUE)
Sample <- args[1]
RNA_Sample <- gsub("MSM","MSN",Sample)

##loading log normalized RNA data including imputed ATAC cells
seRNA <- get(load(paste0(Sample,"_imputed_SCE.RData"))) 

## combined tumor single-cell HVG list
com.hvg <- readLines("Files/topHVGG34new10xnoRPMTSX.txt") #com.hvg

#load post-QC ATAC data, per sample
proj <- loadArchRProject(path = Sample)

DIR = "~/SCENIC/output"
folder_path <- paste0(DIR,"RNA/",Sample,"/")

# Check if the folder already exists
if (!file.exists(folder_path)){
  dir.create(folder_path)
}

setwd(paste0(DIR,"RNA/",Sample,"/"))

#convert RNA SCE to Scanpy AnnData
############ add missing information and filter for ATAC cells ###############

#remove normal cells from sceRNA
remove_cells <- readLines(paste0("normal_cells/",Sample,"_normal_cells.txt"))

print("Number of normal cells")
length(remove_cells)

print("Number of cells before removing Normal cells")
length(colnames(seRNA))

k <- data.frame(X=colnames(seRNA))
k <- k %>% separate(X,c("A","B"))
k2 <- paste0(Sample,"#",k$B,"-1")
colnames(seRNA) <- k2

seRNA <- seRNA[,colnames(seRNA) %ni% remove_cells]

print("Number of cells after removing Normal cells")
length(colnames(seRNA))


#remove cells that are not present in ATAC sample
cellsRNA <- colnames(seRNA)
cellsATAC <- getCellNames(proj)

cellsToKeep <- which(cellsRNA %in% cellsATAC)
print("Number of cells used for SCENIC+:")
length(cellsToKeep)

seRNA <- seRNA[,cellsToKeep]


######## subset RNA genes to contain only HVGs ############
#filter for HVG
#subst the SCE to contain only the top HVG and the cells also present in the ATAC project
seRNA <- seRNA[rownames(seRNA) %in% com.hvg,]

#add the celltype(own annotation) and Combined Cluster information
ATAC <- data.frame(cellname=getCellNames(proj), Clusters_Combined = proj$Clusters_Combined)
RNA <- data.frame(cellname = colnames(seRNA))

clusterinfo <- RNA %>%left_join(ATAC, by= "cellname")

seRNA$Clusters_Combined <- clusterinfo$Clusters_Combined

#save seRNA used for SCENIC+
print("Saving seRNA object used for SCENIC+")
save(seRNA, file = "seRNA_SCENIC.RData")

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
cellMeta<-as.data.frame(seRNA@colData)

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
cellMeta=r.cellMeta
adata.obs=cellMeta

#save adata object
adata.write(os.path.join(projDir, 'adata.h5ad'), compression='gzip')
exit
q()
print(paste0("Saved adata object for sample: ", Sample))
