## This script combines data from all the RNA libraries
suppressPackageStartupMessages({
  library(batchelor)
  library(scran)
  library(scater)
  library(rlist)
})

## merging individual experiments ----
libs <- readLines("Files/MBsamples.txt",header = FALSE)
all.sce <- list()
all.dec <- list()
for(x in libs){
  all.sce[[x]] <- get(load(paste0("files/",x,"scepro.RData")))
  all.dec[[x]] <- get(load(paste0("files/",x,"dec.RData")))
  rm(list=ls(pattern="^SN"))

}

universe <- Reduce(intersect, lapply(all.sce, rownames))
all.sce <- lapply(all.sce, "[", i=universe,)
save(all.sce, file="comdata.RData")

universe2 <- Reduce(intersect, lapply(all.dec, rownames))
all.dec <- lapply(all.dec, "[", i=universe2,)
save(all.dec, file="comdatadec.RData")

# # plotdata ----
sams <- names(all.sce)
cl <- c("nUMIs","nGenes","pct.mt","MALAT1","in_ex_frac","decontX_contamination",
        "SCRscore","sizeFactor","total_counts","S.score","G2M.score","CC.score",
        "CC.per","ind_cluster","Batch","Stage", "Age","Sex",
        "Sample","Region","iUMAP1","iUMAP2")
plot.data <- c()
for(x in sams){
  sce <- get(load(paste0("files/",x,"scepro.RData")))
  rm(list=ls(pattern="^SN"))
  del <- data.frame(colData(sce))
  del2 <- reducedDim(sce, "UMAP")
  del$iUMAP1 <- del2[,1]
  del$iUMAP2 <- del2[,2]

  plot.data <- rbind(plot.data,del[,cl])
  rm(del,del2,sce)
}
save(plot.data, file="plotdataG34MB.RData")

##here I am trying to find 1op HVGs from each and combine them------
# load("comdatadec.RData")
# com.hvg <- c()
# 
# for(x in sams){
#   del <- getTopHVGs(all.dec[[x]])
#   com.hvg <- c(com.hvg,del[1:1500])
# }
# com.hvg <- com.hvg[!duplicated(com.hvg)]
# rm(all.dec)
# 
# load("Files/HGNCsexchr.RData")
# com.hvg <- com.hvg[!(com.hvg %in% SX)]
# write(com.hvg, file="topHVGG34new10xnoRPMTSX.txt")
# q()
#cosine norm----
#sams <- unique(plot.data$Batch)
# com.sce_cnnc <- c()
# for(x in sams){
#   sce <- get(load(paste0("files/",x,"scepro.RData")))
#   rm(list=ls(pattern="^SN"))
#   del <- logcounts(sce)
#   del <- cosineNorm(del[com.hvg,])
#   com.sce_cnnc <- cbind(com.sce_cnnc,del)
#   rm(del,sce)
# }
# save(com.sce_cnnc, com.hvg, file="cosnormG34.RData")
# dim(com.sce_cnnc)
# rm(com.sce_cnnc)
# q()
# # ##------
# load("comdata.RData")
# # #combine counts
# # com.sce <- do.call(noCorrect,list.append(all.sce, assay.type = "counts"))
# # assayNames(com.sce) <- "counts"
# # save(com.sce, file = "combinedcountsG34.RData")
# # rm(com.sce)
# # 
# #MBN
# norm.sce <- do.call(multiBatchNorm,
#                     list.append(all.sce,norm.args=list(use_altexps=FALSE)))
# rm(all.sce)
# 
# com.sce.rnnc <- do.call(noCorrect,norm.sce)
# assayNames(com.sce.rnnc) <- "logcounts"
#mdt <- logcounts(com.sce.rnnc)
#save(mdt, file = "lgcountsG34MBnew10x.RData")

q()