## This script uses output of scrublet filtered doublets and performs final processing
## of the RNA data using scater/scran pipeline
## This script is going to process individual samples for HVG, UMAP, clusters
## and also identify cell cycling score
# 1. Data input----
suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(dplyr)
  library(tidyr)
  
})


##functions ----
##add module score function of Seurat
#obj=logcounts data
# features= list of gene candidates, not array
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
#cell cycle score as percentage of CC gene expression
cc.score <- function(use_genes){
  x <- assay(sce, "counts")
  sf <- colSums(x)
  score <- x[use_genes,]
  sf2 <- colSums(score)
  score <- (sf2/ sf) * 100
    return(score)
}
#argument
args = commandArgs(trailingOnly=TRUE) #sample id
#for MB----
sam <- read.delim("MBsamples.txt", header = TRUE, row.names = 1) ##sample metadata file
load(paste0(args[1],"_postSCR.RData")) ##post scrublet data
##remove doublets from ATAC
load("ATAC_doubelts_to_be_removed_from_RNA.RData") 
Doublets <- df[df$ID==args[1],"Doublets"]
if(length(Doublets)>0){
  keep <- colnames(sce) %in% Doublets
  sce <- sce[,!keep]
}
#scran normalization ----
cl2 <- quickCluster(sce) 
sce <- computeSumFactors(sce, clusters = cl2, min.mean = 0.1 )
#summary(sizeFactors(sce))
sce <- logNormCounts(sce)
sce$total_counts <- colSums(counts(sce))
dd <- as(matrix(0,dim(counts(sce))[1],dim(counts(sce))[2]),"dgCMatrix")
dd[counts(sce)>0]<-1
sce$nGenes <- colSums(dd)
## cell cycle score----
load("Files/SeuratCCgenes.RData")
cc <- read.delim("Files/ccgenesHUM.txt")
cc <- as.character(cc$Gene.ID)
cc <- cc[cc %in% row.names(sce)]
#validated addmodscore by comparing to Seurat's function
del <- addmodscore(logcounts(sce),features = list(s.genes,g2m.genes, cc))
#head(del)
sce$S.score <- del$Cluster1
sce$G2M.score <- del$Cluster2
sce$CC.score <- del$Cluster3
sce$CC.per <- cc.score(cc)

##recalling HVGs and clustering-----

# #removing cell cycle and associated genes
# genes <- row.names(sce)
# mito <- grep(pattern = "^MT-", x = genes,
#             ignore.case = TRUE, value = TRUE)
# rps <- genes[startsWith(genes,"RPS")]
# rpl <- genes[startsWith(genes,"RPL")]
# mrps <- genes[startsWith(genes,"MRPS")]
# mrpl <- genes[startsWith(genes,"MRPL")]
# #ccr <- c(cc,mito, rps,rpl,mrps,mrpl)
# ccr <- c(mito, rps,rpl,mrps,mrpl)
# keep <- genes %in% ccr
# genes <- genes[!keep]
# #save(genes, file= "genesnocc.RData")
# save(genes, file= "genesnoRPMT.RData")

load("genesnoRPMT.RData")
#just getting variation for no cc genes
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
sce$ind_cluster <- paste0(args[1],"_",cluster)

## convert to SCE
##for MB ----
sce$Batch <- args[1]
sce$Group <- as.character(sam[args[1],"Group"])
sce$Subtype <- as.character(sam[args[1],"Subtype"])
sce$Age <- sam[args[1],"Age"]
sce$Sex <- as.character(sam[args[1],"Sex"])
sce$Sample <- as.character(sam[args[1],"Sample"])

rn <- colnames(sce)
rn <- gsub(paste0(args[1],"_"),"",rn)
at <- sam[args[1],"ATAC"]
atc <- paste0(at,"#",rn,"-1")
sce$ATAC_cells <- atc
sce$ATAC_Batch <- at

##remove unwanted meta info
sce$Cluster <- sce$ClusterProb <- sce$score.debris <- sce$Call <-
  sce$Rank <- sce$decontX_clusters <- sce$SCRcall <- NULL
#UMAP
set.seed(408)
del <- uwot::umap(reducedDim(sce,'PCA_sel'), n_neighbors=25, min_dist = 0.5, metric="cosine")
reducedDim(sce,'UMAP') <- del
#cell selection----
row.names(mdfr2) <- paste0(args[1],"_",row.names(mdfr2))
keep.c <- row.names(mdfr2) %in% colnames(sce)
mdfr2$Cell <- "No"
mdfr2[keep.c,"Cell"] <- "Yes"

keep.c <- row.names(mdt.pd) %in% colnames(sce)
mdt.pd$Cell <- "No"
mdt.pd[keep.c,"Cell"] <- "Yes"

##saving ----
assign(args[1],sce)
assign(paste0(args[1],".dec"),sce.dec)
paste0(args[1]) %>% save(list =.,file = paste0(args[1],"scepro.RData"))
paste0(paste0(args[1],".dec")) %>% save(list =.,file = paste0(args[1],"dec.RData"))
#QC plots ------

library(ggsci)
library(ggrepel)
library(gridExtra)
mdt.pd$Cell <- factor(mdt.pd$Cell, levels = c("No","Yes") )
#post selection plot------
#head(mdt.pd)
p1 <- ggplot(mdt.pd%>% arrange(Cell), aes(x=nUMIs, y=nGenes, color= Cell)) +
  geom_point(size=0.2, alpha = 0.5)+
  scale_color_manual(values = c("grey77","red" ))+
  theme_classic()+ xlab("Total Counts") + ylab ("Number Genes") +
  theme(axis.line = element_line(colour = 'black',size=0.5),
        axis.ticks = element_line(colour = 'black',size=0.5),
        axis.text=element_text(size=5),
        axis.title=element_text(size=5,face="bold"),
        legend.position = "None")

p2 <- ggplot(mdt.pd, aes(x=nGenes, y=pct.mt, color= Cell)) +
  geom_point(size=0.2, alpha = 0.5)+
  scale_color_manual(values = c("grey77","red" ))+
  theme_classic()+ xlab("Number Genes") + ylab ("MT%") +
  theme(axis.line = element_line(colour = 'black',size=0.5),
        axis.ticks = element_line(colour = 'black',size=0.5),
        axis.text=element_text(size=5),
        axis.title=element_text(size=5,face="bold"),
        legend.position = "None")

p3 <- ggplot(mdt.pd, aes(x=nGenes, y=MALAT1, color= Cell)) +
  geom_point(size=0.2, alpha = 0.5)+
  scale_color_manual(values = c("grey77","red" ))+
  theme_classic()+ xlab("Number Genes") + ylab ("MALAT1%") +
  theme(axis.line = element_line(colour = 'black',size=0.5),
        axis.ticks = element_line(colour = 'black',size=0.5),
        axis.text=element_text(size=5),
        axis.title=element_text(size=5,face="bold"),
        legend.position = "None")

p4 <- ggplot(mdt.pd, aes(x=pct.mt, y=MALAT1, color= Cell)) +
  geom_point(size=0.2, alpha = 0.5)+
  scale_color_manual(values = c("grey77","red" ))+
  theme_classic()+ xlab("MT%") + ylab ("MALAT1%") +
  theme(axis.line = element_line(colour = 'black',size=0.5),
        axis.ticks = element_line(colour = 'black',size=0.5),
        axis.text=element_text(size=5),
        axis.title=element_text(size=5,face="bold"),
        legend.position = "None")

tiff(paste0(args[1],"_afterem_plot.tiff"),width = 5, height = 5, units ="in", res=300)
grid.arrange(p1,p2,p3,p4, ncol=2)
dev.off()

##fraction plot-----
tiff(paste0(args[1],"_Intronic_frac_Plot.tiff"), width = 4, height = 3, units ="in", res=150)

ggplot(mdfr2, aes(x=Rank, y=frac, color=Cell)) +
  geom_point(size=0.1, alpha=0.5)+
  scale_color_manual(values = c("grey77", "red"))+
  theme_classic()+ labs(y="Fraction of Intronic reads", x="Rank barcodes") +
  theme(axis.line = element_line(colour = 'black',size=0.2),
        axis.ticks = element_line(colour = 'black',size=0.2),
        axis.text = element_text(face="bold",colour = "black",size=5),
        axis.title=element_text(face="bold",colour="black",size = 5),
        panel.grid.major = element_line(colour = "grey97"),
        panel.border =element_rect(fill=NA, colour = 'black',size=0.4),
        legend.position = "None")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+
  scale_x_continuous(expand = c(0, 0))

dev.off()
q()
##UMAP plots ----
plot.data <- data.frame(colData(sce))
del <- data.frame(reducedDim(sce, "UMAP"))
colnames(del) <- c("UMAP1", "UMAP2")
plot.data <- cbind(del,plot.data)
d <- length(unique(plot.data$ind_cluster))
plot.data$ind_cluster <- factor(plot.data$ind_cluster, levels = paste0(args[1],"_",seq(d)))


label.d = plot.data %>% group_by(ind_cluster) %>%
  select(UMAP1, UMAP2) %>% summarize_all(median)
pa <- ggplot(plot.data, aes(x=UMAP1, y=UMAP2, color=ind_cluster))+
  geom_point(size=0.5)+
  scale_color_manual(values =colorRampPalette(pal_simpsons()(14))(d))+
  theme_classic()+
  geom_label_repel(aes(label = ind_cluster),size = 1.5, data = label.d, show.legend = FALSE)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")

pb <- ggplot(plot.data, aes(x=UMAP1, y=UMAP2, color=S.score))+
  geom_point(size=0.5)+
  scale_colour_gradient2(low = "blue", mid= "white", high = "red4")+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")

pc <- ggplot(plot.data, aes(x=UMAP1, y=UMAP2, color=G2M.score))+
  geom_point(size=0.5)+
  scale_colour_gradient2(low = "blue", mid= "white", high = "red4")+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")

pd <- ggplot(plot.data, aes(x=UMAP1, y=UMAP2, color=CC.score))+
  geom_point(size=0.5)+
  scale_colour_gradient2(low = "blue", mid= "white", high = "red4")+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")

tiff(paste0(args[1],"_UMAP.tiff"), width = 8, height = 6, units ="in", res=300)

gridExtra::grid.arrange(pa, arrangeGrob(pb, pc,pd), ncol = 2)

dev.off()

#per cluster QC metrics----
p1 <- ggplot(plot.data, aes(x=ind_cluster, y=decontX_contamination, color=ind_cluster))+
  geom_boxplot(size=0.5)+
  scale_color_manual(values =colorRampPalette(pal_simpsons()(14))(d))+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black',size=0.3),
        axis.ticks = element_line(colour = 'black',size=0.3),
        axis.text.x = element_text(face="bold",colour = "black",size=5, angle = 90, hjust=1, vjust=0.5),
        axis.text.y = element_text(face="bold",colour = "black",size=5),
        axis.title.y=element_text(face="bold",colour="black",size = 5),
        axis.title.x = element_blank(),
        legend.position = "none")


p3 <- ggplot(plot.data, aes(x=ind_cluster, y=SCRscore, color=ind_cluster))+
  geom_boxplot(size=0.5)+
  scale_color_manual(values =colorRampPalette(pal_simpsons()(14))(d))+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black',size=0.3),
        axis.ticks = element_line(colour = 'black',size=0.3),
        axis.text.x = element_text(face="bold",colour = "black",size=5, angle = 90, hjust=1, vjust=0.5),
        axis.text.y = element_text(face="bold",colour = "black",size=5),
        axis.title.y=element_text(face="bold",colour="black",size = 5),
        axis.title.x = element_blank(),
        legend.position = "none")

p4 <- ggplot(plot.data, aes(x=ind_cluster, y=pct.mt, color=ind_cluster))+
  geom_violin(size=0.5)+
  scale_color_manual(values =colorRampPalette(pal_simpsons()(14))(d))+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black',size=0.3),
        axis.ticks = element_line(colour = 'black',size=0.3),
        axis.text.x = element_text(face="bold",colour = "black",size=5, angle = 90, hjust=1, vjust=0.5),
        axis.text.y = element_text(face="bold",colour = "black",size=5),
        axis.title.y=element_text(face="bold",colour="black",size = 5),
        axis.title.x = element_blank(),
        legend.position = "none")

p5 <- ggplot(plot.data, aes(x=ind_cluster, y=in_ex_frac, color=ind_cluster))+
  geom_violin(size=0.5)+
  scale_color_manual(values =colorRampPalette(pal_simpsons()(14))(d))+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black',size=0.3),
        axis.ticks = element_line(colour = 'black',size=0.3),
        axis.text.x = element_text(face="bold",colour = "black",size=5, angle = 90, hjust=1, vjust=0.5),
        axis.text.y = element_text(face="bold",colour = "black",size=5),
        axis.title.y=element_text(face="bold",colour="black",size = 5),
        axis.title.x = element_blank(),
        legend.position = "none")


p6 <- ggplot(plot.data, aes(x=ind_cluster, y=log2(total_counts), color=ind_cluster))+
  geom_violin(size=0.5)+
  scale_color_manual(values =colorRampPalette(pal_simpsons()(14))(d))+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black',size=0.3),
        axis.ticks = element_line(colour = 'black',size=0.3),
        axis.text.x = element_text(face="bold",colour = "black",size=5, angle = 90, hjust=1, vjust=0.5),
        axis.text.y = element_text(face="bold",colour = "black",size=5),
        axis.title.y=element_text(face="bold",colour="black",size = 5),
        axis.title.x = element_blank(),
        legend.position = "none")

tiff(paste0(args[1],"_QC_Plot.tiff"), width = 12, height = 6, units ="in", res=300)

gridExtra::grid.arrange(arrangeGrob(p1,p3),arrangeGrob(p4,p5,p6), nrow=2)

dev.off()

