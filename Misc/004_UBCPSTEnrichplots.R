## plotting gene set enrichment signature along UBC pseudotime
library(ggplot2)
DIR_MB <- "~/MBnewana/"

##auc
load("pdGCUBC_G34MBWGCNA_AUC.RData")
load("pdGCUBC_G34GRN_AUC.RData")
load("pdGCUBC_G34MBNMF_AUC.RData")

##metadata
load("UBC_slingshot_lig5.RData")

##WGCNA
sel <- readLines("module_order.txt")
sel <- gsub("ME","",sel)
##representatives of integrated G34 axes
sel2 <- c("ME16","ME9",
          "ME12","ME2")

##GRN
sel <- read.delim("GRNs.txt")
sel <- sel$GRNs
##selected GRNs
sel2 <- c("NRL","CRX","CREB5","MYC","FOXN4",
          "EOMES","LHX1","ALX1")

#G34 MB NMF
sel2 <- sel <- c("NMF_2","NMF_1")
##sub
sel2 <-sel <- paste0("NMF_",c(3,8,2,7,5,6,4,1))

##########
ms <- scr[row.names(pdUBC),sel]
pdUBC$Rankbin <- as.character(cut_interval(pdUBC$Rank,100))

pd2x <- c()

for(x in sel){
  del <- scale(ms[,x])
  del[del>2] <- 2
  del[del<(-2)] <- (-2)
  del_df <- data.frame(X=pdUBC$Rank,Y=del)
  lmod10 <- loess(Y~X, data=del_df, span=0.2) 
  del <- predict(lmod10)
  pd2x <- cbind(pd2x,del)
}
rm(x)
row.names(pd2x) <- row.names(pdUBC)
rc <- order(pdUBC$Rank)

colnames(pd2x) <- sel
pdUBCx <- pdUBC[rc,]
pd2x <- pd2x[rc,]
cl_lv <- unique(pdUBCx$Rankbin)
msd2 <- c()
for(x in cl_lv){
  keep <- pdUBCx$Rankbin==x
  del <- pd2x[keep,]
  del2 <- apply(del, 2, median)
  msd2 <- rbind(msd2,del2)
  rm(keep,del,del2)
}
rm(x)
row.names(msd2) <- cl_lv
msd2[1:4,]

msd2[msd2>1] <- 1
msd2[msd2<(-1)] <- (-1)

myColor <- colorRampPalette(c("deepskyblue4","beige","magenta3"))(101)
myBreaks <- c(seq(-1, 0, length.out=ceiling(100/2) + 1),
              seq(1/100, 1, length.out=floor(100/2)))

hp <- pheatmap::pheatmap(t(msd2),scale = "row",
                         color=myColor, breaks = myBreaks,
                         cluster_rows = FALSE,cluster_cols = FALSE,
                         show_colnames = FALSE,
                         fontsize_row = 7,
                         border_color = NA,
                         legend = FALSE)

w=length(sel)*0.2+1
tiff("GCUBC_TFaucGRN_stripHM_UBC.tiff", unit="in", width =6 , height = w, res=300) 
print(hp)
dev.off()
