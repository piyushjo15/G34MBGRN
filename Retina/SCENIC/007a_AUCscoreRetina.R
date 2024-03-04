suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(dplyr)
  library(tidyr)
  library(AUCell)
  library(BiocParallel)
  
})
##----
DIR_ATAC <- "Retina/ATAC/"
DIR_RNA <- "Retina/RNA/"
DIR_GRN <- "Retina/ATAC/SCENIC/scenicplus/"
DIR_GRNx <- "SCENIC/GRNana/"

load(paste0(DIR_GRN,"TFselcl_GRN_Retina.RData")) #Retina GEPs
#load(paste0(DIR_GRNx,"TFselcl_GRN_AUC_G34MBv2.RData")) #G34MB GRNs from AUC approach, loading gene sets identifed from combine GRNS

GSEA <- regs
###load  RNA data-----
load(paste0(DIR_RNA,"lgcountsRet.RData"))
load(paste0(DIR_RNA,"plotdataRetANN.RData"))
keep <- plot.data$Annotation=="ND"
table(keep)
table(row.names(plot.data)==colnames(mdt))
plot.data <- plot.data[!keep,]
mdt <- mdt[,!keep]
auc <- AUCell_run(mdt, regs, aucMaxRank=nrow(mdt)*0.10, normAUC = TRUE,
                  BPPARAM=BiocParallel::MulticoreParam(5))
scr <- t(assay(auc)) ## row cell, col GRN
save(plot.data, auc, scr, file =paste0(DIR_GRN,"pdRetR_RetaucGRN_AUC.RData")) ##Retina Retina GRN
#save(plot.data, auc, scr, file =paste0(DIR_GRN,"pdRetR_G34aucGRN_AUC.RData")) ##Retina G34MB GRN

q()
