#preparing each RNA library data for SCENIC
suppressPackageStartupMessages({
  library(scater)
  library(scran)
  library(Matrix)
})
#argument
Samples = commandSamples(trailingOnly=TRUE)

##without normal clusters------
sce <- get(load(paste0(Samples[1],"scepro.RData")))
load("plotdata_G34MB_SVM_Liger.RData")
keep <- plot.data$Nor=="No" ##remove normal cells
plot.data <- plot.data[keep,]

plot.data <- plot.data[plot.data$Batch==Samples[1],]
#for HVG
com.hvg <- readLines("Files/topHVGG34new10xnoRPMTSX.txt")

pc <- read.table("Files/gencdH38p13r37_PC.txt", header = TRUE)
pc <- as.character(pc$Gene.ID)
com.hvg <- com.hvg[com.hvg %in% pc]

exprMat <- round(counts(sce)[com.hvg,row.names(plot.data)])
cellInfo <- data.frame(colData(sce))


library(SCopeLoomR)


SCE <- build_loom(paste0(Samples[1],"_nn.loom"),dgem=exprMat)
SCE <- add_cell_annotation(SCE, cellInfo)
close_loom(SCE)

q()

