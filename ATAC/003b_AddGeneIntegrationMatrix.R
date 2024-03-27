## This script adds predicted gene expression data to each cell from RNA data 
## for the non-multi-omic sample and then obtain joint ATAC+RNA SVD representation
## to obtain UMAP representation and clustering using combined reduced dimensions

#load packages
suppressPackageStartupMessages({
  library(ArchR)
})

#define SampleIDs
Sample <- "SA194"
RNA_Sample <- "SN407"

#setwd
DIR <- "ATAC/"
setwd(paste0(DIR,Sample))
proj <- loadArchRProject(path = Sample)

#load RNA data (dataframe containing all inforamtion + SCE object)
load(paste0("RNA/",RNA_Sample,"scepro.RData"))


########## Add Gene Expresssion Matrix - GeneIntegrationFunction from ArchR #############

sceRNA <- SN407
gr_obj <- rtracklayer::import("Files/gencdH38p13r37CR_genesN.filtered.bed")
gmn <- gr_obj$name
sceRNA <- sceRNA[gmn,]

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneExpressionMatrix",
  reducedDims = "IterativeLSI_qc",
  seRNA = sceRNA,
  addToArrow = TRUE, 
  groupRNA = "ind_cluster",
  nameCell = "predicted_RNA", #store the cell ID from the matched scRNA-seq cell
  nameGroup = "predictedGroup", #store the group ID from the scRNA-seq cell
  nameScore = "predictedScore", #store the cross-platform integration score
  force = TRUE ,
  useImputation = FALSE
)


########## Iterative LSI for RNA data ########################
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


#UMAP and Clustering
proj <- addUMAP(proj, reducedDims = "LSI_RNA", name = "UMAP_RNA_LSI",  minDist = 0.8, force = TRUE)
proj <- addClusters(proj, reducedDims = "LSI_RNA", name = "Clusters_RNA_LSI", resolution = 1, force = TRUE)


############### Joint LSI and UMAP #############################

#get Matrix from proj
mdt <- getMatrixFromProject(proj, "GeneExpressionMatrix")
mdt.M <- assay(mdt)
rownames(mdt.M) <- mdt@elementMetadata@listData$name
#make sce with predcicted GEX
mdt_x <- sceRNA_predicted <- SingleCellExperiment(list(logcounts=mdt.M))
#saving this predicted data for later 
save(mdt_x, file = paste0("Imputation/",Sample,"_imputed.RData"))

#identify HVG
hvg <- getTopHVGs(sceRNA_predicted, n = 1500) ## no need to remove Sex chromosome genes for a sample
mdt_hvg <- mdt.M[rownames(mdt.M) %in% hvg,]

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



saveArchRProject(ArchRProj = proj, load = F)
print("Finished!")
