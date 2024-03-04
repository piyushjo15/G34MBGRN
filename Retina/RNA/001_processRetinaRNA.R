## Processing RNA data for Retina
# 1. Data input----
suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(dplyr)
  library(tidyr)
  library(Seurat)
})
Sample = commandSample(trailingOnly=TRUE)
meta_sams <- read.delim("Files/Ret_metadata.txt",row.names = 1)

mdt <- Read10X(paste0(Sample[1],"/"), gene.column = 1)
##converting ENS id to gene names
ens <- read.delim("Files/gencdH38p13r37CR_genes.txt", header = FALSE)
ens <- ens %>% separate(V1,c("A","B"))
keep <- duplicated(ens$A)
ens <- ens[!keep,]
row.names(ens) <- ens$A

rn <- row.names(mdt) %in% row.names(ens)
rn <- row.names(mdt)[rn]
ens <- ens[rn,]
mdt <- mdt[rn,]
row.names(mdt) <- ens$V2

sce <- SingleCellExperiment(list(counts=mdt))
#scran normalization ----
cl2 <- quickCluster(sce) 
sce <- computeSumFactors(sce, clusters = cl2, min.mean = 0.1 )
#summary(sizeFactors(sce))
sce <- logNormCounts(sce)
sce$total_counts <- colSums(counts(sce))

genes <- row.names(sce)


mito <- grep(pattern = "^MT-", x = genes,
            ignore.case = TRUE, value = TRUE)
rps <- genes[startsWith(genes,"RPS")]
rpl <- genes[startsWith(genes,"RPL")]
mrps <- genes[startsWith(genes,"MRPS")]
mrpl <- genes[startsWith(genes,"MRPL")]
ccr <- c(mito, rps,rpl,mrps,mrpl)
keep <- genes %in% ccr
genes <- genes[!keep]

sce2 <- sce[genes,]
sce.dec <- modelGeneVar(sce2)
top.hvgs <- getTopHVGs(sce.dec, n = 1500)
sce <- runPCA(sce, subset_row = top.hvgs)
rm(sce2)
#selecting number of imp PCA
reducedDim(sce, "PCA_sel") <- reducedDim(sce, "PCA")[,1:15]

#cluster using jaccard-based weights and louvain method as in Seurat
#buildSNNGraph first makes knn graph and then SNN justlike seurat
#heigher k, lower number of clusters
set.seed(126)
g <- buildSNNGraph(sce, k = 20,  use.dimred="PCA_sel", type = "jaccard")
cluster <- igraph::cluster_louvain(g)$membership
sce$ind_cluster <- paste0(Sample[1],"_",cluster)

#UMAP
set.seed(408)
del <- uwot::umap(reducedDim(sce,'PCA_sel'), n_neighbors=25, min_dist = 0.3, metric="cosine")
reducedDim(sce,'UMAP') <- del

##add meta data
sce$Batch <- Sample[1]
sce$Stage <- meta_sams[Sample[1],"Age"]
sce$Sex <- meta_sams[Sample[1],"Sex"]
##saving ----
assign(Sample[1],sce)
assign(paste0(Sample[1],".dec"),sce.dec)
paste0(Sample[1]) %>% save(list =.,file = paste0(Sample[1],"scepro.RData"))
paste0(paste0(Sample[1],".dec")) %>% save(list =.,file = paste0(Sample[1],"dec.RData"))
q()
