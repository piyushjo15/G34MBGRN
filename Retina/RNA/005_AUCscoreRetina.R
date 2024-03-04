## Script to calculate AUC scores in Retna
suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(dplyr)
  library(tidyr)
  library(AUCell)
  library(BiocParallel)
  
})
##part 1 ..module score ----
load("TFselcl_GRN_AUC_G34MBv2.RData") #G34MB GRNs 

###load  RNA data-----
load("lgcountsRet.RData")
load("plotdataRetANN.RData")
keep <- plot.data$Annotation=="ND" ## remove low quality cells
table(keep)
table(row.names(plot.data)==colnames(mdt))
plot.data <- plot.data[!keep,]
mdt <- mdt[,!keep]
genes <- row.names(mdt)

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


auc <- AUCell_run(mdt, GSEA, aucMaxRank=nrow(mdt)*0.10, normAUC = TRUE,
                  BPPARAM=BiocParallel::MulticoreParam(5))
scr <- t(assay(auc)) ## row cell, col GRN
save(plot.data, auc, scr, file =paste0(DIR_GRN,"pdRetR_G34aucGRN_AUC.RData")) ##Retina G34MB GRN

q()
#### Enrichment heatmap #######
##tf enrichment by ind_cluster per batch
zscore_c <- function(x){
  sdx <- sd(x)
  mx <- mean(x)
  ot <- boxplot.stats(x)$out
  q <- boxplot.stats(x)$stats
  keep <- (x < q[1])
  x[keep] <- q[1]
  keep <- (x > q[5])
  x[keep] <- q[5]
  x2 <- x-mx
  z <- x2/sdx
  return(z)
}

ms2 <- apply(scr, 2, zscore_c)
msd2 <- c()
cl_lv <- c("Pro_early","Pro_late1","Pro_late2","Pro_late_II1","Pro_late_II2",
           "RGC_e1","RGC_e2","RGC","AC/HC_pro","AC_e","AC_1","AC_2",
           "HC_e1","HC_e2","HC",
           "Cones_e1","Cones_e2","Cones_l",
           "Rod_e","Rods_l","Rods_II",
           "BiP_e1","BiP_e2","BiP_on1","BiP_on2","BiP_on3","BiP_off1","BiP_off2","BiP_off3",
           "Astro","MG")

for(x in cl_lv){
  keep1 <- plot.data$Annotationlv2==x
  msx <- ms2[keep1,]
  del <- colMeans(msx)
  msd2 <- rbind(msd2,del)
  rm(keep1,msx,del)
}
rm(x)
rownames(msd2) <- cl_lv
msd2[1:4,1:4]
sel <- read.delim("GRNs.txt")
sel <- sel$GRNs

myColor <- colorRampPalette(c("deepskyblue4","beige","magenta3"))(101)
myBreaks <- c(seq(-2, 0, length.out=ceiling(100/2) + 1),
              seq(2/100, 2, length.out=floor(100/2)))
####-----

hp <- pheatmap::pheatmap(t(msd2[,sel]), 
                         cluster_cols = FALSE, 
                         cluster_rows = FALSE,
                         #show_rownames = FALSE,
                         legend = FALSE,
                         color=myColor, breaks=myBreaks,
                         border_color = NA,
                         fontsize = 8)
tiff("Ret_G34MBTFGRN_enr.tiff", width = 6, height = 12, units = "in", res = 300)
print(hp)
dev.off()
