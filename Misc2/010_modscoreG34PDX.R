##This scripts calculates enrichment score of gene-sets in combined tumor data
suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(AUCell)
  library(BiocParallel)
})
##----

load(paste0(DIR_GRN,"TFselcl_GRN_AUC_G34MBv2.RData")) 

###load  RNA data-----
load("lgcountsG34MBPDX.RData")
load("plotdata_G34MB_PDX.RData")
table(row.names(plot.data)==colnames(mdt))
###AUCell-----
ms <- AUCell_run(mdt, GSEA, aucMaxRank=nrow(mdt)*0.10, normAUC = TRUE,
                 BPPARAM=BiocParallel::MulticoreParam(9))
scr <- t(assay(ms))
save(plot.data, ms, file =paste0(DIR_GRN,"pdHDMB03_G34GRN_AUCv2.RData")) ##G34MB G34MB GRN

q()
