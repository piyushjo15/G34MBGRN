#preparing data for SCENIC
suppressPackageStartupMessages({
  library(scater)
  library(scran)
  library(Matrix)
})
#argument
args = commandArgs(trailingOnly=TRUE)
##without normal clusters------
sce <- get(load(paste0(args[1],"scepro.RData")))
colnames(sce) <- paste0(args[1],"_",colnames(sce))
load("plotdataRetANN.RData")
keep <- plot.data$Annotation=="ND" ##remove abnormal cells
plot.data <- plot.data[!keep,]

plot.data <- plot.data[plot.data$Batch==args[1],]
#for HVG
load("topHVGRetRiboMTSX.RData")
com.hvg <- topHVG[1:5000]

pc <- read.table("Files/gencdH38p13r37_PC.txt", header = TRUE)
pc <- as.character(pc$Gene.ID)
com.hvg <- com.hvg[com.hvg %in% pc]

exprMat <- round(counts(sce)[com.hvg,row.names(plot.data)])
cellInfo <- plot.data


library(SCopeLoomR)


SCE <- build_loom(paste0(args[1],"_ret.loom"),dgem=exprMat)
SCE <- add_cell_annotation(SCE, cellInfo)
close_loom(SCE)

q()
