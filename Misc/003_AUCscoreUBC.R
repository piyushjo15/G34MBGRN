## AUC score of gene-sets in cerebellar UBC lineage
suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(dplyr)
  library(tidyr)
  library(AUCell)
  library(BiocParallel)
  
})

load("referencedata.RData")
## CBsce.rnnc is the combined SCE data and plot.data is the metadata

##subsetting to GC/UBC lineage, with GC-defined
sel <- c("progenitor","GCP/UBCP","GCP","GC_diff_1","GC_diff_2",
         "UBC_diff", "UBC_defined")
keep <- plot.data$dev_state %in% sel
plot.data <- plot.data[keep,]
CBsce.rnnc <- CBsce.rnnc[,row.names(plot.data)]
common.genes <- readLine() ## top HVG in GC/UBC lineage 

genes <- common.genes[1:10000]
load("TFselcl_GRN_AUC_G34MBv2.RData") #G34MB GRNs from AUC app
#load("GS_WGCNA.RData") ##MB WGCNA
#load("G34MB_NMF_metagene.RData") ##G34MB NMF metagene
#GSEA <- G34MBSub
# ##-----
#some filtering
GS <- GSEA
sets <- names(GS)
GSEA <- list()
for(x in sets){
  del <- GS[[x]]
  keep <- del %in% genes
  del <- del[keep]
  GSEA[[x]] <- del
  rm(del)
}
rm(x)


########AUC############ does not look good############

auc <- AUCell_run(mdt, GSEA,
                  aucMaxRank=nrow(mdt)*0.10, normAUC = TRUE,
                  BPPARAM=BiocParallel::MulticoreParam(8))
scr <- t(assay(auc)) ## row cell, col GRN
save(pdUBC, auc,scr, file ="pdGCUBC_G34GRN_AUC.RData") ## G34Mb TF-GRN
#save(pdUBC,auc,scr, file = "pdGCUBC_G34MBWGCNA_AUC.RData") ##WGCNA programs 
#save(pdUBC,auc,scr, file = "pdGCUBC_G34MBNMF_AUC.RData") ##G34MB NMF metagene
q()