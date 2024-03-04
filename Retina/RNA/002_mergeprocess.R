## Merging all retina RNA library
suppressPackageStartupMessages({
  library(batchelor)
  library(scran)
  library(scater)
  library(rlist)
})

# ## merging individual experiments ----
libs <- readLines("Retinalibs.txt")
all.sce <- list()
all.dec <- list()
plot.data <- c()
for(x in libs){
  sce <- get(load(paste0(x,"scepro.RData")))
  colnames(sce) <- paste0(x,"_",colnames(sce))
  all.sce[[x]] <- sce
  all.dec[[x]] <- get(load(paste0(x,"dec.RData")))
  del <- data.frame(colData(sce))
  del2 <- reducedDim(sce, "UMAP")
  del$iUMAP1 <- del2[,1]
  del$iUMAP2 <- del2[,2]
  plot.data <- rbind(plot.data,del)
  
  rm(sce,del,del2)

}
rm()
save(plot.data, file="plotdataRet.RData")

universe <- Reduce(intersect, lapply(all.sce, rownames))
all.sce <- lapply(all.sce, "[", i=universe,)
save(all.sce, file="comdata.RData")

universe2 <- Reduce(intersect, lapply(all.dec, rownames))
all.dec <- lapply(all.dec, "[", i=universe2,)
save(all.dec, file="comdatadec.RData")

com.dec <- do.call(combineVar,all.dec)
rm(all.dec)
topHVG <- getTopHVGs(com.dec)
load("Files/HGNCsexchr.RData")
topHVG <- topHVG[!(topHVG %in% SX)]
save(topHVG, file="topHVGRetRiboMTSX.RData")
com.sce_cnnc <- c()
for(x in libs){
  sce <- all.sce[[x]]
  del <- logcounts(sce)
  del <- cosineNorm(del[topHVG[1:4000],]) 
  com.sce_cnnc <- cbind(com.sce_cnnc,del)
  rm(del,sce)
}
save(com.sce_cnnc, topHVG, file="cosnorRet.RData")
rm(com.sce_cnnc)
com.sce <- do.call(noCorrect,list.append(all.sce, assay.type = "counts"))
assayNames(com.sce) <- "counts"
save(com.sce, file = "combinedcountsRet.RData")
rm(com.sce)


norm.sce <- do.call(multiBatchNorm,
                    list.append(all.sce,norm.args=list(use_altexps=FALSE)))
rm(all.sce)

com.sce.rnnc <- do.call(noCorrect,norm.sce)
assayNames(com.sce.rnnc) <- "logcounts"

mdt <- logcounts(com.sce.rnnc)
save(mdt, file = "lgcountsRet.RData")

q()