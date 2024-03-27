## This script add gene expression data to each cell from its multi-omic RNA 
## counterpart. If the cells passed QC in ATAC but not in RNA, then we impute
## log-counts data for those ATAC-cells. We also obtain joint RNA+ATAC SVD reduced 
## dimension for  UMAP representing and clustering

library(reticulate)

Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")

#load packages
suppressPackageStartupMessages({
  library(ArchR)
  library(scran)
  library(scater)
  library(Matrix)
  library(batchelor)
})

set.seed(456)

#setwd
addArchRGenome("hg38")
Sample <- commandArgs(TRUE) ##These are ids from MBataclibs.txt except SA194
RNA_Sample <- gsub("MSA","MSN",Sample)
DIR="ATAC/"
setwd( paste0(DIR,Sample,"/"))

##load RNA data
sceRNA <- get(load(paste0(RNA_Sample,"scepro.RData")))
plot.data_RNA <- data.frame(colData(sceRNA))

## load ATAC metadata generated post Preprocessing. In some samples, few clusters
## that exhibited low fragment number and TSS enrichment were removed. This was 
## done manually
load(paste0(DIR,Sample,"/",Sample,"_barcode_stats.RData"))

#convert Clusters_qc to characters as previously this was converted into factors
plot.data$Clusters_qc <- as.character(plot.data$Clusters_qc)
#
plot.data$RNA <- FALSE
cells_rna <- plot.data_RNA$ATAC_cells ## cells passing QC in RNA analysis
keep <- row.names(plot.data) %in% cells_rna #define the ATAC cells that are also present in RNA data
print("Number of overlapping RNA-ATAC cells: ")
print(sum(keep))

plot.data[keep,"RNA"] <- TRUE #define in ATAC plot.data which ATAC cells have RNA cells or not
plot.data$Impute <- !(plot.data$RNA) #make a new column in plot.data with the opposite entry of "RNA" - cells that require imputation

## count number of cells need to be imputed per ATAC cluster
cl_imp <- data.frame(table(plot.data$Clusters_qc,plot.data$RNA))
cl_imp <- reshape(cl_imp, idvar = "Var1", timevar = "Var2", direction = "wide") # separate column for FALSE and TRUE
colnames(cl_imp) <- c("Clusters_qc","No","Yes")
row.names(cl_imp) <- cl_imp$Clusters_qc
cl_imp$Clusters_qc <- NULL
cl_imp$Total <- cl_imp$No+cl_imp$Yes

#convert to fraction
cl_imp$No <- (cl_imp$No/cl_imp$Total)*100
cl_imp$Impute <- "Yes"
cl_imp$Remove <- "No"
cl_imp[cl_imp$No>50,"Impute"] <- "No" #don't impute these cluster cells, minimum 50% have RNA data
cl_imp[cl_imp$Yes<100,"Remove"] <- "Yes" ##remove all the cells from this cluster, low number of reference cell
##
print(cl_imp)

#remove cells
cl_to_rem <- row.names(cl_imp)[cl_imp$Remove=="Yes"]
if(length(cl_to_rem)>0){plot.data <- plot.data[!(plot.data$Clusters_qc %in%cl_to_rem),]}
cl_to_not_impute <- row.names(cl_imp)[cl_imp$Impute=="No"]
if(length(cl_to_not_impute)>0){plot.data[(plot.data$Clusters_qc %in%cl_to_not_impute),"Impute"] <- FALSE} #set Impute column to FALSE

print("Number of cells after filtering: ")
print(length(plot.data$Clusters_qc))

#save plot.data
save(plot.data,file=paste0(DIR,Sample,"/",Sample,"_barcode_stats_Imputation.RData")) #plot.data


print(paste0("Finished Filtering",Sample))


################ ATAC: load the data ########################
proj <- loadArchRProject(path = Sample)
#load the plot.data for the information which cells require imputation
plot.data_Imputation <- get(load(paste0(DIR,Sample,"/",Sample,"_barcode_stats_Imputation.RData"))) 

top.hvgs <- getTopHVGs(sceRNA, n = 1500) ## no need to remove Sex chromosome genes for a sample

############# prepare ATAC part for Imputation #################
#load LSI matrix from ATAC data
LSI_qc <- getReducedDims(proj,reducedDims = "IterativeLSI_qc")
LSI_qc[1:4,1:4]
cell_atac <- row.names(LSI_qc)

############# prepare RNA part for Imputation ###########################

mdt <- logcounts(sceRNA)
cells_rna <- plot.data_RNA$ATAC_cells
colnames(mdt) <- cells_rna
plot.data_RNA$Imputed <- "No"
plot.data_RNA$ATAC <- "No"

##find intersecting cells
cells_com <- intersect(cell_atac, cells_rna)
cells_4imp <- rownames(plot.data_Imputation[plot.data_Imputation$Impute == TRUE,]) ## cells with only good atac data

##define which cells have to be removed from ATAC
remove <- rownames(plot.data_Imputation[plot.data_Imputation$RNA == plot.data_Imputation$Impute,]) #remove cells that have no RNA data and are not imputed
remove2 <- cell_atac[cell_atac %ni% rownames(plot.data_Imputation)] #remove the cells that have been removed during the filtering script

print("Number of RNA cells: ")
length(colnames(sceRNA))

print("Number of RNA-ATAC intersecting cells: ")
length(cells_com)

print("Number of cells to impute: ")
length(cells_4imp)

print("Number of cells to remove from ATAC: ")
length(c(remove,remove2))

#remove
keep <- which(getCellNames(proj) %ni% c(remove,remove2))
subsetArchRProject(proj,
                   cells = getCellNames(proj)[keep],
                   outputDirectory = Sample,
                   force=TRUE)

#define lsi_train and lsi_test
lsi_train <- LSI_qc[cells_com,] #we know RNA for this data
lsi_test <- LSI_qc[cells_4imp,] # we need to find RNA for these cells

mdt_train <- mdt[,cells_com] ##the RNA for train data
plot.data_RNA[plot.data_RNA$ATAC_cells %in% cells_com,"ATAC"] <- "Yes"


################# Imputation ##############################
#using KNN to identify nearest neighbors
## python code ----
repl_python()

import numpy as np
import sklearn.neighbors as skn

neigh = skn.NearestNeighbors(n_neighbors=5) #reduced from 100 to 25 since no huge difference could be seen, later to 5 for better distribution
neigh.fit(r.lsi_train)
dis, nd = neigh.kneighbors(r.lsi_test, return_distance=True) ##get distance for weighted KNN

exit

##converting into R objects
kneigh <- py$nd ##index of neighbors
kneigh <- kneigh + 1 ##add 1 to index, python vs R differencce
kneigh_dist <- py$dis ##index of neighbors
row.names(kneigh) <- row.names(kneigh_dist) <- cells_4imp
colnames(kneigh) <- colnames(kneigh_dist) <-paste0("N_",seq(dim(kneigh)[2]))

head(kneigh)
head(kneigh_dist)

##converting distance matrix for weight matrix
max_d <- apply(kneigh_dist, 1, max)
kneigh_dist_i <- 1/kneigh_dist #inversing the distance to contribution, farther has lower weight
sum_dis <- rowSums(kneigh_dist_i)
kneigh_w <-  kneigh_dist_i/sum_dis##normalizing  
head(kneigh_w)

##finding psuedocounts from neighbors
mdt_quer <- as(matrix(0,dim(mdt_train)[1],dim(lsi_test)[1]),"dgCMatrix")
row.names(mdt_quer) <- row.names(mdt_train)
colnames(mdt_quer) <- cells_4imp
mdt_quer[1:4,1:4]

for(x in cells_4imp){
  n_x <- kneigh[x,]
  #weighted
  n_w <- kneigh_w[x,]
  mdt_r <-t(t(mdt_train[,n_x])*n_w)##obtain data from neighbors, multiply each sample by weights and transpose
  mdt_quer[,x] <- rowSums(mdt_r) ##sum all the expression
  rm(n_x,n_w, mdt_r)

}
rm(x)

mdt_x <- cbind(mdt,mdt_quer)
mdt_comb <- cbind(mdt_train,mdt_quer)
mdt_comb[1:4,1:4]

###create meta data for all RNA cells + imputed
pd <- plot.data_RNA[1:length(cells_4imp),]
k <- data.frame(X=cells_4imp)
k <- k %>% separate(X,c("A","B"))
k2 <- paste0(RNA_Sample,"_",k$B)

head(k2)

row.names(pd) <- k2
pd[,c(1:14,21)] <- NA #set it to NA
pd$Imputed <- pd$ATAC <- "Yes"
pd$ATAC_cells <- cells_4imp
head(pd)

plot.data_x <- rbind(plot.data_RNA,pd) #plot.data = all RNA cells, pd = imputed cells
save(mdt_x,mdt,mdt_quer, plot.data, plot.data_x, file = paste0("Imputation/",Sample,"_imputed.RData"))

#save imputed RNA as SCE object
sceRNA_imputed <- SingleCellExperiment(list(logcounts=mdt_x))

##adding colData info to new SCE object
#first check order of the cells
plot.data_cells <- row.names(plot.data_x)
k <- data.frame(X=plot.data_cells)
k <- k %>% separate(X,c("A","B"))
k2 <- paste0(Sample,"#",k$B,"-1")
plot.data_cells <- k2
table(plot.data_cells==colnames(sceRNA_imputed))
#add colData (plot.data_x)
colData(sceRNA_imputed) <- DataFrame(plot.data_x)
save(sceRNA_imputed, file = paste0("Imputation/",Sample,"_imputed_SCE.RData"))

###create meta data for only those RNA cells where ATAC data is present + imputed = plot.data_comb
plot.data_comb <- plot.data_x[plot.data_x$ATAC == "Yes",]
rownames(plot.data_comb) <- plot.data_comb$ATAC_cells

#add Imputed information to proj - First checking the order of the cells to be the same in plot.data_comb and mdt_comb with ArchR proj!!!!
colData_ATAC <- getCellColData(proj)
cellsATAC <- rownames(colData_ATAC)
cellsImp <- rownames(plot.data_comb)

# Get the order 
order <- match(cellsATAC, cellsImp)
# Reorder vector 1
plot.data_comb <- plot.data_comb[order,]
cellnames <- rownames(plot.data_comb)

cellsRNA <- colnames(mdt_comb)
cellsATAC <- getCellNames(proj)

#order the colnames of mdt_hvg to be in the same order as ATAC!
# Get the order 
order <- match(cellsATAC, cellsRNA)
# Reorder vector 1
mdt_comb <- mdt_comb[,order]
cellnames <- colnames(mdt_comb)

print("Compare the order of barcodes between ATAC cells and cellnames of RNA mdt_comb matrix:")
is_same_order <- identical(cellnames, cellsATAC)
# Print the result
if (is_same_order) {
  print("The order of barcodes is the same.")
} else {
  print("The order of barcodes is different.") #the order is the same!
}

########### add imputed information to the proj + SCE cluster information ##################
proj$Imputed <- plot.data_comb$Imputed
proj$SCE_Cluster <- plot.data_comb$ind_cluster
proj$nGenes <- plot.data_comb$nGenes
proj$nUMIs <- plot.data_comb$nUMIs
proj$total_counts <- plot.data_comb$total_counts


######################### Add Gene ExpressionMatrix ######################################
#create Summarized Experiment from the imputed data
#mdt_x = normalized GEX matrix
gr_obj <- rtracklayer::import("Files/gencdH38p13r37CR_genesN.filtered.bed")

#check if the rownames of plot.data is in the same order as the colnames of mdt_comb - it should be the same, just checking..
cell_mdt <- colnames(mdt_comb)
cell_plot <- rownames(plot.data_comb)

is_same_order <- identical(cell_mdt,cell_plot)

print("checking order of barcodes in mdt_comb and plot.data_comb for generating SE....")
# Print the result
if (is_same_order) {
  print("The order of barcodes is the same.")
} else {
  print("The order of barcodes is different.") #the order is the same!
}

se <- SummarizedExperiment(assays = list(logcounts=mdt_comb), colData = plot.data_comb,
                           rowRanges = gr_obj)

#convert cell names to the same as proj

proj <- addGeneExpressionMatrix(
  input = proj,
  seRNA = se,
  chromSizes = getChromSizes(proj),
  excludeChr = c("chrM", "chrY")
)


################ Multiome  ###################
############# take SCRAN normalized matrix ##############

#define hvg for the sample (already defined above)
hvg <- top.hvgs

#subst the mdt_comb to contain only the top HVG
mdt_hvg <- mdt_comb[rownames(mdt_comb) %in% hvg,]

cellsRNA <- colnames(mdt_hvg)
cellsATAC <- getCellNames(proj)
#order the colnames of mdt_hvg to be in the same order as ATAC!
# Get the order 
order <- match(cellsATAC, cellsRNA)
# Reorder vector 1
mdt_hvg <- mdt_hvg[,order]
cellnames <- colnames(mdt_hvg)

print("Compare the order of barcodes between ATAC and RNA matrix for Joint LSI:")
is_same_order <- identical(cellnames, cellsATAC)
# Print the result
if (is_same_order) {
  print("The order of barcodes is the same.")
} else {
  print("The order of barcodes is different.") #the order is the same!
}

# Perform SVD on the mdt_hvg
library(SingleCellExperiment)
library(irlba)

dim <- 30
sce_svd_hvg <- irlba(mdt_hvg, nv = dim) #used as input for svd for LSI_RNA_scran! (top HVG : dim)

# Extract SVD results
svd_matrix_u <- sce_svd_hvg$u 
svd_matrix_v <- sce_svd_hvg$v #matrix of right singular vectors
svd_values <- sce_svd_hvg$d #vector of singular values

dim(svd_matrix_u) # dimension is top HVG : 30 (nv)
dim(svd_matrix_v) # dimension is cells : 30 (nv)

d <- vector <- paste0("LSI", seq(1, 30))
colnames(svd_matrix_v) <- d
rownames(svd_matrix_v) <- colnames(mdt_hvg)

#make a list for the IterativeLSI_RNA input
#LSI-RNA: run IterativeLSI on the GeneExpressionMatrix to have the list as template

proj <- addIterativeLSI(
  ArchRProj = proj,
  clusterParams = list(
    resolution = 0.2,
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix",
  depthCol = "Gex_nUMI",
  varFeatures = 1500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA",
  force = TRUE
)

ArchR <- proj@reducedDims@listData$LSI_RNA #use this as template

proj@reducedDims@listData$IterativeLSI_RNA <- ArchR
proj@reducedDims@listData$IterativeLSI_RNA$matSVD <- svd_matrix_v #enter the cells:dim matrix
proj@reducedDims@listData$IterativeLSI_RNA$svd <- sce_svd_hvg

#Add UMAP and Clusters
proj <- addUMAP(proj, reducedDims = "IterativeLSI_RNA", name = "UMAP_RNA",  minDist = 0.8, force = TRUE)
proj <- addClusters(proj, reducedDims = "IterativeLSI_RNA", name = "Clusters_RNA", resolution = 1, force = TRUE)

################   combine IterativeLSI_RNA + IterativeLSI_qc ######################
proj <- addCombinedDims(proj, reducedDims = c("IterativeLSI_qc", "IterativeLSI_RNA"), name =  "LSI_Combined")

#UMAP
proj <- addUMAP(proj, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)

#Clustering
proj <- addClusters(proj, reducedDims = "LSI_Combined", name = "Clusters_Combined", resolution = 1, force = TRUE)

#save ATAC proj and new barcode-stats
saveArchRProject(ArchRProj = proj, load = F)

plot.data <- getCellColData(proj)
save(plot.data,file=paste0(DIR,"Samples_Ensembl/",Sample,"/",Sample,"_barcode_stats_afterMulti.RData"))
print("saved proj after adding LSI_Combined")
## generate a table for 1 column = cell barcode + 1 column = cluster IDs (clusters_combined) + 1column = imputed, yes or no?

table <- data.frame(barcode = as.character(getCellNames(proj)), Clusters_Combined = as.character(proj$Clusters_Combined), Imputed = as.character(proj$Imputed))
write.table(table,paste0(Sample,"_cellname_cluster_df.txt"), sep = "\t", row.names = F, col.names = T)

print(paste0("saved after Multiome Script for Sample: ", Sample))



