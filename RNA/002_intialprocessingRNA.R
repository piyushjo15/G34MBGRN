## This script uses GeneFull (intronic) and Gene (exonic) counts from STAR aligner to generate a combined gene matrix
## The combined gene matrix is used by diem to identify cell vs debris.
## The output of diem is used by SoupX and decontX to remove debris.
## The final output is saved 
## SoupX is done on fixed counts, which includes preMRNA + exons counts of overlapping/undercounted genes
# 1. Data input----
suppressPackageStartupMessages({
  library(diem)
  library(ggplot2)
  library(gridExtra)
  library(Matrix)
  library(scater)
  library(scran)
})

#argument
args = commandArgs(trailingOnly=TRUE) ##Sample id

## on server----
#Load the premrna counts
DIR="~/path/to/rawcounts/"
counts.in <- read_10x(paste0(DIR,args[1],"/GeneFull"))
dim(counts.in)
# load the exon part
counts.ex <- read_10x(paste0(DIR,args[1],"/Gene"))
dim(counts.ex)

## intronic read fraction plot------
tot <- colSums(counts.in)
exonic <- colSums(counts.ex)
#reads not present in exonic reads are intronic?
intronic <- tot-exonic  
frac <- intronic/tot
mdfr <- data.frame(cbind(frac,tot))

#filtering Nan
mdfr <- mdfr[!(is.na(frac)),]
#sorting by tot reads
mdfr <- mdfr[order(-mdfr$tot),]
#rank
mdfr$Rank <- seq(1:dim(mdfr)[1])
##subsetting 20k cells
mdfr2 <- mdfr[1:20000,]
rm(mdfr, intronic, exonic, tot, frac)

#check colnames
table(colnames(counts.in) ==colnames(counts.ex))
#check rownames
table(row.names(counts.in) == row.names(counts.ex))
#check example
head(counts.ex)[,1:10]
genes <- row.names(counts.in)
#finding genes that were counted less with intronic marix
mean.in <- rowMeans(counts.in)
mean.ex <- rowMeans(counts.ex)
d <- (mean.ex > mean.in)
table(d)
#replacing above genes' data in intronic matrix from exonic matrix
counts.in2 <- counts.in[!d,]
dim(counts.in2)
counts.ex2 <- counts.ex[d,]
dim(counts.ex2)

counts.in2 <- rbind(counts.in2,counts.ex2)
dim(counts.in2)
counts.in2 <- counts.in2[genes,]
rm(genes)

#check colnames and rownames
if(all(colnames(counts.in) ==colnames(counts.ex)) & all(row.names(counts.in) == row.names(counts.ex))){
  
  #head(counts.ex)[,1:10]
  genes <- read.delim("~/HB/ann/gencdH38p13r37CR_genesN.txt", header = FALSE) ##ENS id
  genes <- unique(genes$V1)
  
  #finding genes that were counted less with intronic marix
  mean.in <- rowMeans(counts.in)
  mean.ex <- rowMeans(counts.ex)
  d <- (mean.ex > mean.in)
  print(table(d))
  #replacing above genes' data in intronic matrix from exonic matrix
  counts.in2 <- counts.in[!d,]
  dim(counts.in2)
  counts.ex2 <- counts.ex[d,]
  dim(counts.ex2)
  
  counts.in2 <- rbind(counts.in2,counts.ex2)
  dim(counts.in2)
  counts.in2 <- counts.in2[genes,]
  rm(genes)
}

## 1.1 Creating an SCE object------

mdt <- create_SCE(counts.in2, name=args[1])
dim(mdt)
#class(mdt)
### remove all count matrices
rm(counts.ex2, counts.ex, counts.in)
#processing using diem

## 1.2 Plotting quality metrics
mt_genes <- grep(pattern = "^MT-", x = rownames(mdt@gene_data),
                 ignore.case = TRUE, value = TRUE)
mdt <- get_gene_pct(x = mdt, genes = mt_genes, name = "pct.mt")
genes <- grep(pattern = "^MALAT1$", x = rownames(mdt@gene_data),
              ignore.case = TRUE, value = TRUE)
mdt <- get_gene_pct(x = mdt, genes = genes, name = "MALAT1")


# head(drop_data)
# summary(drop_data)
#
# 2. Running DIEM-----

# The following outlines the steps involved in `diem`:
#
# 1. Speciyfing the test set and the debris set.
# 2. Running PCA on test set
# 2. Running k-means on the PCs to initialize the clusters.
# 4. Running a EM to estimate the parameters of the Dirichlet-multinomial mixture model.
# 5. Classifying droplets based on their likelihood.

## 2.1 Specifying test and debris droplets
mdt <- set_debris_test_set(mdt, top_n = 15000,
                           min_counts = 200,
                           min_genes = 200)

length(mdt@test_set)
length(mdt@bg_set)

## 2.2 Running PCA on the test set and intialization-----

mdt <- filter_genes(mdt, cpm_thresh = 10)
genes <- gene_data(mdt)
#summary(genes)
##default pca parameters, genes=2000, pcs=30
mdt <- get_pcs(mdt, threads = 8)


## 3.3 Initializing clusters
#`k_init = 30` says to initialize k-means with 30 centers,
#`min_size_init = 20` says to take clusters with less than 20 droplets
#and assign them to the next closest cluster.
#`nstart_init = 50`  specify to run k-means 50 times and pick the best run.
#`iter.max_init` says each run can have at most 100 iterations.
# We run multiple starts because k-means may converge to a local optimum.
# The run with the lowest total within sum of
# squares is selected.
mdt <- init(mdt,
            k_init = 30,
            nstart_init = 50,
            min_size_init = 50,
            seedn = 123,
            threads = 8)

## 3.4 Running EM using DM like before-----
mdt <- run_em(mdt, threads=8)
mdt <- assign_clusters(mdt)
#less max_genes make debris score lower, use default settings
mdt <- estimate_dbr_score(mdt)

#de_genes <- debris_genes(mdt)
#head(de_genes)
mdt <- call_targets(mdt,thresh_score = .5, min_genes = 200)##
mdt.pd <- droplet_data(mdt)
colnames(mdt.pd)[1:2] <- c("nUMIs", "nGenes")
mdt.pd <- mdt.pd[order(-mdt.pd$nUMIs),]
mdt.pd$Rank <- seq(dim(mdt.pd)[1])
mdt.pd$nUMIs <- mdt.pd$nUMIs+1 ##to account for zero

# save ----
save(mdt, mdfr2,mdt.pd, file=paste0("diemv4",args[1],".RData"))
#load(paste0("diemv4",args[1],".RData"))
cells <- row.names(mdt.pd[mdt.pd$Call=="Clean",]) ##filter for "Clean" cells
sce <- SingleCellExperiment(list(counts=mdt@counts[,cells])) ##conver to SCE object
#table(row.names(drop_data)==colnames(sce))
del <- mdt.pd[colnames(sce),]
colData(sce) <- DataFrame(del)

mdfr3 <- mdfr2[colnames(sce),]
sce$in_ex_frac <- mdfr3$frac
##add sample id to cell names
row.names(mdt.pd) <- paste0(args[1],"_",row.names(mdt.pd))
colnames(sce) <- paste0(args[1],"_",colnames(sce))

#head(colData(sce))
#initial filtering based on pct.mt and feature counts,
# already filtered based on min_gene =200 in diem run
x <- mean(sce$pct.mt)
z <- x+ 1*mad(sce$pct.mt)
if (z <2) z =2
##based on previous data 6600 genes are in 98% quantile
keep <- (sce$pct.mt < z) & (sce$nGenes < 6600)
sce <- sce[,keep]
keep <- sce$in_ex_frac >0.25 #intronic fraction more than 25%
sce <- sce[,keep]

rm(x,keep,z, mdt, mdfr3)

#removing soup using SoupX----
suppressPackageStartupMessages({
  library(SoupX)
  library(Seurat)
  library(dplyr)
})

#to select important debris genes
#Converting into Seurat object for running SoupX
sce.S <- CreateSeuratObject(counts = counts(sce))
sce.S <- SCTransform(sce.S, verbose = FALSE)
sce.S <- RunPCA(sce.S)
sce.S <- RunUMAP(sce.S, dims = 1:30)
sce.S <- FindNeighbors(sce.S, dims = 1:30)
sce.S <- FindClusters(sce.S, resolution = 1)

filt.matrix <- counts(sce)
soup.channel  <- SoupChannel(counts.in2, filt.matrix) ##counts.in2 is the fixed intronic gene expression matrix, line 64
rm(counts.in2)
meta    <- sce.S@meta.data
umap    <- sce.S@reductions$umap@cell.embeddings
soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)
#head(meta)
soup.channel  <- autoEstCont(soup.channel, tfidfMin = 0.8, soupQuantile = 0.85, forceAccept = TRUE)

soup.data <- soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ]

soup.data$symbol <- row.names(soup.data)

tiff(paste0(args[1],"_SoupX_geneplot.tiff"),width = 6, height = 4, units ="in", res=150)

ggplot(filter(soup.data, rank(-est) < 100) , aes(x=rank(-est), y=est, label=symbol)) +
  geom_point(size=0.1, alpha=0.3) +
  geom_text(check_overlap = T) +
  theme_bw()

dev.off()


soupx_corrected  <- adjustCounts(soup.channel, roundToInt = T)
assay(sce,"unfiltcounts") <- counts(sce)
counts(sce) <- soupx_corrected

##running decontX
library(celda)
sce <- decontX(sce, z = sce.S$seurat_clusters)
assay(sce,"rawcounts") <- counts(sce)
counts(sce) <-  decontXcounts(sce)
decontXcounts(sce) <- NULL

save(sce,mdfr2,mdt.pd, file = paste0(args[1],"_postSoupXDX.RData")) ##This will be used to remove doublets

q()

