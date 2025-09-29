## This script uses snATAC-seq data to obtain inferCNV profile per sample
## the input is cellranger output for atac: 
## atac_fragments.tsv.gz  atac_fragments.tsv.gz.tbi  filtered_feature_bc_matrix.h5
## and cell annotation with reference cells annotated

## 1. obtain annotation file 
# ## for  run on axial identities using ATAC cells---------
sam <- "MB129"
## Obtain normal cells for a sample, and add it to tumor cells which are annotated
## for axis or cell-state
load("plotdata_G34MB_integrated.RData")## integrated tumor data with normal
load("plotdata_G34MBaucKNN.RData")## intergrated tumor data without normal

## obtain normal cells
pd <- plot.data[plot.data$Batch==sam,] #subset for batch

pd <- pd[!(pd$Nor=="ND"),]
pd_nor <- row.names(pd)[pd$Nor=="Yes"]
ann_nor <- data.frame(Cluster=rep("Normal",length(pd_nor)),Class=rep("Normal",length(pd_nor)))
row.names(ann_nor) <- pd_nor
head(ann_nor)

## obtain tumor cells and annotation
pdx <- plot.data.com[plot.data.com$Batch==sam,] #subset for batch
pdx <- pdx[pdx$ATAC=="Yes",]
pdx <- pdx[!(pdx$Axis=="MYC"),] ## remove MYC cells which are rare for MB129,MB292 and MB26

ann <- pdx[,c("Annotation","Axis")]
colnames(ann) <- c("Cluster","Class")
head(ann)
table(row.names(ann) %in% row.names(ann_nor)) ##mutually exclusive
ann <- rbind(ann,ann_nor)
rn <- row.names(ann)
rn <- gsub(paste0(sam,"_"),"",rn)
rn <- paste0(rn,"-1")
row.names(ann) <- rn
head(ann)
table(ann$Cluster,ann$Class)

write.table(ann, file=paste0("CNV/",sam,"_ann_ax.txt"), sep="\t",quote = FALSE)##atac v3

## 2. run atac infer CNV---------
suppressPackageStartupMessages({
  library(atacInferCnv)
  library(infercnv)
})

mainDir = "CNV/"
setwd(mainDir)
inDir = paste0(mainDir,sam)
sId = paste0(sam,"_sub")
sAnn = paste0(mainDir,sam ,"_ann_ax.txt" ) ## annotation from ATAC cells
resPath=paste0(mainDir,sId,"_result/")
prepareAtacInferCnvInput(inDir,sAnn,resPath,targColumn = "Class",ctrlGrp = "Normal")

runAtacInferCnv(resPath, numClusters = 1, analysis_mode="samples") # numClusters=3 for MB26, 2 for MB129 and 1 for MB292
# numClusters were obtained from empirical selection with clusters ranging from 1-4

