## Deconvolution analysis of bulk data, dataset specific
## to avoid batch effect I am dividing samples by dataset

suppressPackageStartupMessages({
  library(scater)
  library(scran)
  library(dplyr)
  library(tidyr)
  library(BayesPrism)
  library(DESeq2)
})
DIR_MB <- "RNA/"
DIR_GRN <- "SCENIC/GRNana/"
##processing dataset for HVG----
load("MDT_DS2.RData") ## or ICGC or CR
keep <- MDT$Subgroup_M %in% c("G3","G4")
table(keep)
MDT <- MDT[keep,]

mt <- tvg.MDT[startsWith(tvg.MDT,"MT-")]
rps <- tvg.MDT[startsWith(tvg.MDT, "RPS")]
mrps <- tvg.MDT[startsWith(tvg.MDT, "MRPS")]
rpl <- tvg.MDT[startsWith(tvg.MDT, "RPL")]
mrpl <- tvg.MDT[startsWith(tvg.MDT, "MRPL")]
#removed human mitochondria list
rmg <- c(rps,rpl,mt,mrpl,mrps)

keep <- tvg.MDT %in% rmg
tvg.MDT <- tvg.MDT[!keep]
load("Files/HGNCsexchr.RData")
tvg.MDT <- tvg.MDT[!tvg.MDT %in% SX]

##reference----
com.hvg <- readLines("Files/topHVGG34new10xnoRPMTSX_1500.txt")

##getting normal cells----------
## this is integrated tumor data with non-neuronal normal cells
load(paste0(DIR_MB,"plotdata_G34MB_integrated.RData")) 
plot.data <- plot.data[plot.data$Nor=="Yes",]
table(plot.data$SVM_prilabKD)

### Merging certain normal cell-types into a class
##astroglia 
rn_astro <- row.names(plot.data)[plot.data$SVM_prilabKD%in% c("astrocytes","progenitor")]
astro <- data.frame(Annotation=rep("Astroglia",length(rn_astro)),KNN_cl=rep("Astroglia",length(rn_astro)))
row.names(astro) <- rn_astro
##oligodendro 
rn_oligo <- row.names(plot.data)[plot.data$SVM_prilabKD %in% c("oligo_progenitor","oligodendrocyte")]
oligo <- data.frame(Annotation=rep("Oligodendrocytes",length(rn_oligo)),KNN_cl=rep("Oligodendrocytes",length(rn_oligo)))
row.names(oligo) <- rn_oligo
##immune
rn_imm <- row.names(plot.data)[plot.data$SVM_prilabKD=="immune"]
imm <- data.frame(Annotation=rep("Immune",length(rn_imm)),KNN_cl=rep("Immune",length(rn_imm)))
row.names(imm) <- rn_imm
##endo 
rn_endo <- row.names(plot.data)[plot.data$SVM_prilabKD%in% c("mural/endoth","meningeal")]
endo <- data.frame(Annotation=rep("Endothelial",length(rn_endo)),KNN_cl=rep("Endothelial",length(rn_endo)))
row.names(endo) <- rn_endo

nor <- rbind(astro,oligo,imm,endo)
nor$Axis <- "Normal"
rm(plot.data)
##tumor
load(paste0(DIR_GRN,"plotdata_G34MBaucKNN.RData")) ### THis is tumor data annotated as per 4 axes
plot.data.com <- plot.data.com[!(is.na(plot.data.com$ind_cluster)),]
plot.data.com <- plot.data.com[,c("Annotation","KNN_cl","Axis")]
plot.data.com <- rbind(plot.data.com,nor)


###load tumor raw-counts
load(paste0(DIR_MB,"combinedcountsG34new10x.RData"))
com.sce <- com.sce[,row.names(plot.data.com)]
# #v2
com.hvg <- intersect(com.hvg, tvg.MDT[1:7500])

#need to be matrix
com.sce <- as.matrix(counts(com.sce[com.hvg,row.names(plot.data.com)]))
mdt.MDT <- mdt.MDT[com.hvg,row.names(MDT)]

#need to be transposed
com.sce <- t(com.sce)
mdt.MDT <- t(mdt.MDT)

#v2
ann1 <- plot.data.com$Annotation
ann2 <- plot.data.com$KNN_cl

##run prism----
myPrism <- new.prism(
  reference=com.sce, 
  mixture=mdt.MDT,
  input.type="count.matrix", 
  cell.type.labels = ann1, 
  cell.state.labels = ann2,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1)

theta <- get.fraction (bp=bp.res,
                       which.theta="final",
                       state.or.type="type")
save(theta, MDT,file = "MDT_thetaGT_G34MBANN.RData")
q()
