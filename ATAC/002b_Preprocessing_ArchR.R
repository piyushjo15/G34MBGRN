## This script processes arrow file for doublet identification, clustering and
## intial QC
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyr)
  library(dplyr)
  library(BSgenome.Hsapiens.UCSC.hg38) 
})

set.seed(456)
addArchRThreads(threads =6) 
quant_cut <- function(x){
  min_q <-quantile(x,0.05)
  max_q <-quantile(x,0.95)
  x[x<min_q] <- min_q
  x[x>max_q] <- max_q
  return(x)
}

#Next, I load the reference genome
addArchRGenome("hg38")
Sample <- commandArgs(TRUE) ## This is MBataclibs.txt

DIR="ATAC"

if (dir.exists(paste0(DIR,Sample))==F){
  dir.create(paste0(DIR,Sample))
}

setwd( paste0(DIR,Sample,"/"))

# Removing doublets------
print(paste0("processing sample ",Sample," for ArchR pre-processing step!"))
ArrowFiles <- paste0("ATAC/ArrowFiles/", Sample,".arrow")

if (dir.exists(paste0("doublets"))==F){
  dir.create(paste0("doublets"))
}

print("calculating doublet scores .....")
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "LSI", 
  LSIMethod = 1,
  dimsToUse = 1:50,
  nTrials = 10,
  scaleDims = F,
  threads = 1,
  LSIParams = list(seed = 1,
                   varFeatures = 50000,
                   excludeChr = c("chrX", "chrY", "chrMT")),
  outDir = "doublets"
)
print("creating ArchR project post doublet scores .....")

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = Sample, 
  copyArrows = TRUE 
)

##loading meta data
met_sams <- read.delim("Files/MBsamples.txt",row.names = 2)

proj$Group <- met_sams[Sample,"Group"]
proj$Subtype <- met_sams[Sample,"Subtype"]
proj$Age <- met_sams[Sample,"Age"]
proj$Sex <- met_sams[Sample,"Sex"]
proj$Sample <- met_sams[Sample,"Sample"]
proj$RNA_ID <- meta_sams[Sample,"RNA_ID"]
#number of cells post intial QC
print("Number of cells before doublet filtering:")
nCells(proj)
plot.data <- data.frame(getCellColData(proj))
plot.data$barcode <- row.names(plot.data)


# 1. Filtering putative doublets (pass 1).
# Here we want to be relatively lenient, so we select a ratio of 1 (i.e. filtering top 5% of barcodes for a sample of 5,000 cells).
# We are using Doublet enrichment as a cutoff.

proj <- filterDoublets(ArchRProj = proj,filterRatio = 1)


#save a list storing barcodes of doublets, for this I compare the cells from proj to proj2 and see which cells are missing in proj2
print("Number of cells post doublet filtering:")
nCells(proj)
##saving doublet cells
doublets <- plot.data$barcode[!(plot.data$barcode %in% getCellNames(proj))]
write(doublets, paste0("doublets/",Sample,"/",Sample,"_doublets.txt"))

#2. Filter fragments with high fragment number + cells with high Doublet enrichment
x <- mean(plot.data$nFrags)
z <- x+ 2*mad(plot.data$nFrags)
keep1 <- plot.data$nFrags > z

#remove doublet enrichment
keep2 <- plot.data$DoubletEnrichment > 4
put_doublets <- plot.data$barcode[keep1 | keep2]
length(put_doublets)

#total putative doublets that will be removed
cells_rem <- unique(c(put_doublets,doublets))
write(cells_rem,  paste0("doublets/",Sample,"_rem_cells.txt"))
keep <- plot.data$barcode %in% cells_rem
plot.data <- plot.data[!keep,]
##subsetting object to filtered cells
proj <- proj[plot.data$barcode,]

print("Number of cells post 2nd QC Filtering:")
nCells(proj)


#3. Remove doublest identified in RNA sample

RNA_Sample <- gsub("MSM", "MSN", Sample)
RNA_doublets_file <- paste0(RNA_Sample, "_doublets.txt")

# Check if the RNA_doublets.txt file exists
if (file.exists(RNA_doublets_file)) {
  # Read the RNA_doublets.txt file
  setwd("/omics/odcf/analysis/OE0290_projects/mb_10x/ATACana/RNAdoublets/")
  RNA_doublets <- readLines(RNA_doublets_file)
  
  # Set working directory for the current sample
  setwd(paste0(DIR, "Samples_Ensembl/", Sample, "/"))
  
  # Change to the same barcode format as ATAC data
  RNA_doublets_remove <- paste0(Sample, "#", RNA_doublets, "-1")
  print(paste0("Remove following barcodes as RNA doublets from proj: ", RNA_doublets_remove))
  
  # Check in ATAC data and remove doublets
  cellsToKeep <- which(getCellNames(proj) %ni% RNA_doublets_remove)
  print(paste0("number of cells after removing all doublets:", length(cellsToKeep)))
  proj <- proj[cellsToKeep,]
  
  # Also remove RNA doublets from plot.data
  keep <- plot.data$barcode[plot.data$barcode %ni% RNA_doublets_remove]
  plot.data <- plot.data[keep,]
  
  # Updating these barcodes in our log file
  print("Number of cells post RNA doublet filtering:")
  nCells(proj)
} else {
  print(paste0("RNA_doublets.txt not found for sample ", Sample, ". Continuing without removing RNA doublets."))
}


#Now we can calculate our first iterative LSI. We are now extending our dimensions to 100 and the number of variable features to 100K.

proj <- addIterativeLSI(ArchRProj = proj,
                        useMatrix = "TileMatrix",
                        name = "IterativeLSI_qc",
                        iterations=5,
                        clusterParams = list(
                          resolution = c(0.1, 0.2, 0.4, 0.8),
                          sampleCells = 20000, 
                          n.start = 10
                        ),
                        varFeatures = 100000,
                        dimsToUse = 1:100,
                        totalFeatures = 500000,
                        seed = 1,
                        LSIMethod = 1,
                        scaleDims = FALSE,
                        corCutOff = 0.75,
                        excludeChr = c("chrX", "chrY", "chrMT"),
                        binarize = T,
                        force = T)

##iidentifying clusters to perform second QC pass with high resolution=2
print("Adding clusters:")

proj <- addClusters(input = proj,
                    name = "Clusters_qc",
                    reducedDims = "IterativeLSI_qc",
                    method = "Seurat",
                    force = T,
                    resolution=2,
                    corCutOff = 0.75,
                    scaleDims = FALSE,
                    seed = 1)


# UMAP calculation

proj <- addUMAP(ArchRProj = proj,
                name = "UMAP_qc",
                reducedDims = "IterativeLSI_qc",
                minDist = 0.2,
                metric = "cosine",
                nNeighbors = 25,
                force = T,
                seed = 1,
                scaleDims = F,
                corCutOff = 0.75)


saveArchRProject(proj, outputDirectory = Sample)

save(plot.data, file =paste0(Sample,"_barcode_stats.RData"))

print(paste0("Finished Preprocessing for Sample: ", Sample))
