##This script utilzies RNA data to generate input loom file for pyscenic run
## Normal cells are remvoed and gene expression matrix is subsetted to combined tumor HVG
#preparing data for SCENIC
suppressPackageStartupMessages({
  library(scater)
  library(scran)
  library(Matrix)
})
#argument
Sample = commandArgs(trailingOnly=TRUE)
DIR_RNA <- "RNA/"
##without normal clusters------
#loading RNA data without imputation
sce <- get(load(paste0(DIR_RNA,args[1],"scepro.RData")))
## combined tumor single-cell HVG list
com.hvg <- readLines("Files/topHVGG34new10xnoRPMTSX.txt") #com.hvg
#remove normal cells from sceRNA
remove_cells <- readLines(paste0(DIR_RNA,"normal_cells/",Sample,"_normal_cells.txt"))

keep <- colnames(sce) %in% remove_cells
sce <- sce[,!keep] #subset
plot.data <- data.frame(colData(sce))

pc <- read.table("Files/gencdH38p13r37_PC.txt", header = TRUE)
pc <- as.character(pc$Gene.ID)
com.hvg <- com.hvg[com.hvg %in% pc]

exprMat <- round(counts(sce)[com.hvg,row.names(plot.data)])
cellInfo <- data.frame(colData(sce))


library(SCopeLoomR)


SCE <- build_loom(paste0(args[1],"_nn.loom"),dgem=exprMat)
SCE <- add_cell_annotation(SCE, cellInfo)
close_loom(SCE)

q()
