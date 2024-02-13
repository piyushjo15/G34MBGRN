##This scripts calculates enrichment score of gene-sets in combined tumor data
suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(AUCell)
  library(BiocParallel)
})
##part 1 ..module score ----
addmodscore <- function(obj, features, name = "Cluster",
                        pool = row.names(obj), nbin=24,
                        ctrl = 100, k =FALSE, seed = 123){
  # Find how many gene lists were provided. In this case just one.
  cluster.length <- length(x = features)
  # For all genes, get the average expression across all cells (named vector)
  data.avg <- Matrix::rowMeans(x = obj[pool, ])
  # Order genes from lowest average expression to highest average expression
  data.avg <- data.avg[order(data.avg)]

  # Use ggplot2's cut_number function to make n groups with (approximately) equal numbers of observations.
  # The 'rnorm(n = length(data.avg))/1e+30' part adds a tiny bit of noise to the data, presumably to break ties.
  data.cut <- ggplot2::cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30,
                                  n = nbin,
                                  labels = FALSE,
                                  right = FALSE)
  # Set the names of the cuts as the gene names
  names(x = data.cut) <- names(x = data.avg)
  # Create an empty list the same length as the number of input gene sets.
  # This will contain the names of the control genes
  ctrl.use <- vector(mode = "list", length = cluster.length)

  # For each of the input gene lists:
  for (i in 1:cluster.length) {
    # Get the gene names from the input gene set as a character vector
    features.use <- features[[i]]

    # Loop through the provided genes (1:num_genes) and for each gene,
    # find ctrl (default=100) genes from the same expression bin (by looking in data.cut):
    for (j in 1:length(x = features.use)) {
      # Within this loop, 'data.cut[features.use[j]]' gives us the expression bin number.
      # We then sample `ctrl` genes from that bin without replacement and add the gene names to ctrl.use.
      ctrl.use[[i]] <- c(ctrl.use[[i]],
                         names(x = sample(x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
                                          size = ctrl,
                                          replace = FALSE)))
    }
  }
  # Remove any repeated gene names - even though we set replace=FALSE when we sampled genes from the same expression bin,
  # there may be more than two genes in our input gene list that fall in the same expression bin, so we can end up sampling
  # the same gene more than once.
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)


  ## Get control gene scores

  # Create an empty matrix with dimensions;
  # number of rows equal to the number of gene sets (just one here)
  # number of columns equal to number of cells in input Seurat obj
  ctrl.scores <- matrix(data = numeric(length = 1L),
                        nrow = length(x = ctrl.use),
                        ncol = ncol(x = obj))

  # Loop through each provided gene set and add to the empty matrix the mean expression of the control genes in each cell
  for (i in 1:length(ctrl.use)) {
    # Get control gene names as a vector
    features.use <- ctrl.use[[i]]
    # For each cell, calculate the mean expression of *all* of the control genes
    ctrl.scores[i, ] <- Matrix::colMeans(x = obj[features.use,])
  }

  ## Get scores for input gene sets
  # Similar to the above, create an empty matrix
  features.scores <- matrix(data = numeric(length = 1L),
                            nrow = cluster.length,
                            ncol = ncol(x = obj))

  # Loop through input gene sets and calculate the mean expression of these genes for each cell
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- obj[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }

  # Subtract the control scores from the feature scores -
  # the idea is that if there is no enrichment of the genes in the geneset in a cell,
  # then the result of this subtraction should be ~ 0
  features.scores.use <- features.scores - ctrl.scores

  # Name the result the "name" variable + whatever the position the geneset was in the input list, e.g. "Cluster1"
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)

  # Change the matrix from wide to long
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))

  # Give the rows of the matrix, the names of the cells
  rownames(x = features.scores.use) <- colnames(x = obj)

  return(features.scores.use)

}

##----
DIR_GRN <- "SCENIC/GRNana/"
DIR_BRS <- "BulkRS/"
DIR_MB <- "RNA/"

load(paste0(DIR_GRN,"TFselcl_GRN_AUC_G34MB.RData")) #G34MB GRNs from AUC approach, loading gene sets identifed from combine GRNS
#load(paste0(DIR_MB,"GS_WGCNA.RData")) ##MB WGCNA

###load  RNA data-----
load("ICGC_DS2.RData")
load("CR_DS2.RData")
load("MDT_DS2.RData")
rm(list=ls(pattern = "^dds*"))
rm(mdt.ICGC, mdt.CR, mdt.MDT, meta)
genes <- rownames(mdt.vsd.ICGC)
mdt <- cbind(mdt.vsd.ICGC,mdt.vsd.CR,mdt.vsd.MDT)
meta <- rbind(ICGC,CR,MDT)

##fix GSEA
GSEAx <- list()
sams <- names(GSEA)
for(x in sams){
  del <- GSEA[[x]]
  del <- intersect(del,genes)
  GSEAx[[x]] <- del
  rm(del)
}
GSEA <- GSEAx


ms <- addmodscore(mdt,features = GSEA)
colnames(ms) <- names(GSEA)
#save(meta, ms, file =paste0(DIR_BRS,"pdMB_G34aucGRN_MSv2.RData")) ##G34MB G34MB GRN
save(meta, ms, file =paste0(DIR_BRS,"pdMB_WGCNA_MS2.RData")) ##G34MB WGCNA 

rm(ms)
#q()
###AUCell-----
ms <- AUCell_run(mdt, GSEA, aucMaxRank=nrow(mdt)*0.10, normAUC = TRUE,
                 BPPARAM=BiocParallel::MulticoreParam(9))
scr <- t(assay(ms))
#save(meta, ms, file =paste0(DIR_BRS,"pdMB_G34aucGRN_AUCv2.RData")) ##G34MB G34MB GRN
save(meta, ms,scr, file =paste0(DIR_BRS,"pdMB_WGCNA_AUC.RData")) ##G34MB WGCNA 

q()