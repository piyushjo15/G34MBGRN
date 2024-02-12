##This script calculates enrichment scores per GRN using AUCell, and then identifies
## GRNs that are marker for tumor clusters
suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(dplyr)
  library(tidyr)
  library(AUCell)
  library(BiocParallel)
})

##Directories 
Sample = commandSample(trailingOnly=TRUE)
DIR_GRNana <- "SCENIC/GRNana/"
DIR_ATAC <- "ATAC/"
DIR_RNA <- "RNA/"
DIR_GRN <- "SCENIC/scenicplus/"
msn <- gsub("MSM","MSN",Sample[1])
if(Sample=="SA194"){msn <- "SN407"}
# 
##load expression data----
### load RNA data, without imputation
load(paste0(DIR_RNA,Sample,"scepr0.RData")) 

#remove normal cells from sceRNA
remove_cells <- readLines(paste0(DIR_RNA,"normal_cells/",Sample,"_normal_cells.txt"))
keep <- colnames(sce) %in% remove_cells
sce <- sce[,!keep]
pd1 <- data.frame(colData(sce))
#load GRNs detected for this specific sample
load(paste0(DIR_GRN,Sample[1],"/output/metadata/",msn,"_gmt.RData")) 

cl_lv <- unique(pd1$ind_cluster) ##RNA clusters

##AUC score----
auc <- AUCell_run(mdt, regs, aucMaxRank=nrow(mdt)*0.10, normAUC = TRUE,
                 BPPARAM=BiocParallel::MulticoreParam(5))
scr <- t(assay(auc)) ## row cell, col GRN
save(auc, scr, pd1, file = paste0(DIR_GRNana,msn,"_AUC_GRN.RData")) 

##using pair-wise comaprison from scar::findMarkers, identify marker
## GRNS enriched per cluster. using wilcox to rank enrichment 
marks <- findMarkers(t(scr), groups=pd1$ind_cluster, direction="up", ##row gene
                     test.type="wilcox",pval.type="some")
ups <-c()
all_reg <-c()
##identify top 5 TF-GRNs per cluster
for(x in cl_lv){
  del <- data.frame(marks[[x]])
  del <- row.names(del)[1:3]
  ups <- rbind(ups,del)
  all_reg <- c(all_reg,del)
  rm(del)
}
rm(x)
row.names(ups) <- cl_lv
all_reg <- unique(all_reg)
save(ups, all_reg, file = paste0(DIR_GRNana,msn,"_top3GRN_AUC.RData")) #_AUC for AUC


q()
