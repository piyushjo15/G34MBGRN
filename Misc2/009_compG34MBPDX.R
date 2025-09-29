## This script merges snRNA-seq data from HDMB03 CTRL and PAX6 PDX

suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(tidyverse)
  library(rlist)
  library(batchelor)
})

setwd("RNA/")


## merging PDX data. This is obtained using post-diem, SoupX, scruble, scran/scater processing
load("procdata/HDMB03_CTRLscepro.RData")
load("procdata/HDMB03_PAX6scepro.RData")


load("procdata/HDMB03_CTRLdec.RData")
load("procdata/HDMB03_PAX6dec.RData")

all.sce <- list(HDMB03_CTRL=HDMB03_CTRL,HDMB03_PAX6=HDMB03_PAX6,)

all.dec <- list(HDMB03_CTRL=HDMB03_CTRL.dec,HDMB03_PAX6=HDMB03_PAX6.dec)
sams <- names(all.sce)
plot.data <- c()
for(x in sams){
  del <- data.frame(colData(all.sce[[x]]))
  del2 <- reducedDim(all.sce[[x]], "UMAP")
  del$iUMAP1 <- del2[,1]
  del$iUMAP2 <- del2[,2]
  plot.data <- rbind(plot.data,del)
  rm(del,del2)
}
rm(x)
save(plot.data, file="plotdata_G34MB_PDX.RData")
##
com.hvg <- c()

sams <- names(all.dec)
for(x in sams){
  del <- getTopHVGs(all.dec[[x]])
  com.hvg <- c(com.hvg,del[1:1500])
}
com.hvg <- com.hvg[!duplicated(com.hvg)]

rm(all.dec)

load("HGNCsexchr.RData")
com.hvg <- com.hvg[!(com.hvg %in% SX)]

save(com.hvg, file="topHVG_G34MBPDXnoccSX_1500.RData")

com.sce_cnnc <- c()
for(x in names(all.sce)){
  del <- logcounts(all.sce[[x]])
  del <- cosineNorm(del[com.hvg,])
  com.sce_cnnc <- cbind(com.sce_cnnc,del)
  rm(del)
}
save(com.sce_cnnc, com.hvg, file="cosnormG34MBPDX.RData")

com.sce <- do.call(noCorrect,list.append(all.sce, assay.type = "counts"))
assayNames(com.sce) <- "counts"
save(com.sce, file = "combinedcountsG34MBPDX.RData")
# rm(com.sce)

#MBN
norm.sce <- do.call(multiBatchNorm,
                    list.append(all.sce,norm.args=list(use_altexps=FALSE)))
rm(all.sce)

com.sce.rnnc <- do.call(noCorrect,norm.sce)
assayNames(com.sce.rnnc) <- "logcounts"

mdt <- logcounts(com.sce.rnnc)
save(mdt, file = "lgcountsG34MBPDX.RData")

