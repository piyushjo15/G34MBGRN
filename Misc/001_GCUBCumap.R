## This script extracts counts for GC/UBC lineage cells and LIGER corrects them
suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(batchelor)
  library(rlist)
  library(Matrix)
  library(uwot)
  library(BiocParallel)
  library(liger)
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

##removing low batch samples
del <- data.frame(table(plot.data$batch))
keep <- del$Freq <20

b <- as.character(del[keep,"Var1"])
keep <- plot.data$batch %in% b
plot.data <- plot.data[!keep,]
CBsce.rnnc <- CBsce.rnnc[,row.names(plot.data)]


common.genes <- common.genes[1:1000] 
CBsce.rnnc<- CBsce.rnnc[common.genes,]
# ##liger ------
Tum <- list()
Tum2 <- list()
sams <- unique(plot.data$batch)
for (x in sams) {
  keep <- plot.data$batch==x
  sce <- CBsce.rnnc[,keep]
  del <- logcounts(sce)
  Tum[[x]] <- round(del)
  del <- cosineNorm(del)
  Tum2[[x]] <- as.matrix(t(del))
  rm(sce,del)
}

## Liger
com <- createLiger(Tum)
com@scale.data <- Tum2

com <- optimizeALS(com, k = 15, max.iters=1000) 
del <- com@H
rd_als <- c()
for (x in sams) {
  rd_als <- rbind(rd_als,del[[x]])
}

set.seed(198)
del_umap <- umap(rd_als,n_neighbors = 10, min_dist = 0.1)
plot.data$lUMAP1 <- del_umap[,1]
plot.data$lUMAP2 <- del_umap[,2]

set.seed(777)
mnn.nmf <- reducedMNN(rd_als, batch=plot.data$batch, BPPARAM = MulticoreParam(6))
rd_als_cor <- mnn.nmf$corrected

set.seed(111)
del_umap <- umap(rd_als_cor,n_neighbors =25, min_dist = 0.2, metric="cosine")
plot.data$lcUMAP1 <- del_umap[,1]
plot.data$lcUMAP2 <- del_umap[,2]

save(rd_als, rd_als_cor, plot.data, file= "GCUBCumapLiger5.RData")
q()
