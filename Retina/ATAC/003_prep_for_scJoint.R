## This script integrates Retina ATAC data with Retina RNA data using scJoint

suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(Matrix)
})
DIR_RNA <- "Retina/RNA/"
DIR_ATAC <- "Retina/ATAC/"
DIR_scJoint <- "Retina/ATAC/scJoint/data/"

######################reference data################################
#load reference data

load(paste0(DIR_ATAC,"plotdataRetANN.RData"))
plot.data.ref <- plot.data
rm(plot.data)
#remove ND
keep <- plot.data.ref$Annotation=="ND"
table(keep)
plot.data.ref <- plot.data.ref[!keep,]
load(paste0(DIR_RNA, "lgcountsRet.RData"))
load(paste0(DIR_RNA, "topHVGRetRiboMTSX.RData"))
table(row.names(plot.data.ref) %in% colnames(mdt))
mdt_ref <- mdt[,row.names(plot.data.ref)] ##no duplicate genes
rm(mdt)
######################target data################################
suppressPackageStartupMessages({
  library(batchelor)
  library(scran)
  library(scater)
  library(rlist)
})
# ## merging individual experiments GSM----
libs <- readLines("Files/Retinalibs.txt")
mdt.com <- c()
mdt.com_counts <- c()

for(x in libs){
  load(paste0(DIR_ATAC,x,"/",x,"_GSM.RData"))
  mdt.com_counts <- cbind(mdt.com_counts,mdt)
  sce <- SingleCellExperiment(list(counts=mdt))
  cl2 <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters = cl2, min.mean = 0.1 )
  sce <- logNormCounts(sce)
  mdt.com <- cbind(mdt.com,logcounts(sce))
  rm(mdt,sce,cl2)
}
rm(x)

mdt <- mdt.com
save(mdt, file = paste0(DIR_ATAC,"lgcountsRet_GSM.RData"))
save(mdt.com_counts, file =paste0(DIR_ATAC, "countsRet_GSM.RData"))
#load  data
#load(paste0(DIR_ATAC,"countsRet_GSM.RData"))

##obtained from merging metadata of each Retina ATAC object
load(paste0(DIR_ATAC,"plotdataRet_ATAC.RData")) 
plot.data.tar <- plot.data
rm(plot.data)
##combined plot data
plot.data <- plyr::rbind.fill(plot.data.ref,plot.data.tar)
plot.data$Platform <- "RNA"
keep <- is.na(plot.data$Batch)
table(keep)
plot.data[keep,"Platform"] <- "ATAC"

at <- plot.data[keep,"Sample"]
plot.data[keep,"Batch"] <- paste0(at,"_ATAC")
at <- plot.data[keep,"predictedGroup"]
plot.data[keep,"Annotationlv2"] <- at
row.names(plot.data) <- c(row.names(plot.data.ref),row.names(plot.data.tar))
save(plot.data, file = "com_RNA_ATAC_pd.RData")
rn <- duplicated(row.names(mdt.com_counts))
mdt_tar <- mdt.com_counts[!rn,]
rm(mdt.com_counts)
topHVG <- intersect(topHVG, intersect(row.names(mdt_ref), row.names(mdt_tar)))

source("data_to_h5.R") ##scJoint script
write_h5_scJoint(exprs_list = list(rna = mdt_ref[topHVG[1:5000],],
                                   atac = mdt_tar[topHVG[1:5000],]), 
                 h5file_list = c(paste0(DIR_scJoint,"Retina_rna.h5"), 
                                 paste0(DIR_scJoint,"Retina_atac.h5")))
write_csv_scJoint(cellType_list =  list(rna = plot.data.ref$Annotationlv2,
                                        atac = plot.data.tar$predictedGroup),
                  csv_list = c(paste0(DIR_scJoint,"Retina_cellType_rna.csv"),
                               paste0(DIR_scJoint,"Retina_cellType_atac.csv")))
q()