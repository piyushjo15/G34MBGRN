
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyr)
  library(dplyr)
  library(ggsci)
  library(ggrepel)
  library(tibble)
  library(hexbin)
  library(BSgenome.Hsapiens.UCSC.hg38) 
})

set.seed(456)
addArchRThreads(threads =6) 
Sample = commandSample(trailingOnly=TRUE)
quant_cut <- function(x){
  min_q <-quantile(x,0.05)
  max_q <-quantile(x,0.95)
  x[x<min_q] <- min_q
  x[x>max_q] <- max_q
  return(x)
}
#Next, I load the reference genome
addArchRGenome("hg38")

# Removing doublets------
# For this we will be using simulated doublets and identifying their neighbours in the LSI space.
# Throughout this project, we are using the TF-(logIDF) LSI method as introduced by Cusanovich et al.
# We don't scale LSI components (i.e. each is proportional to its variance) but remove components with more than 0.75 correlation with sequencing depth.
# Here, we will be using 50 dimensions, 10 iterations and label the 10 NNs to each doublet.
##already created Arrow file before
print(paste0("processing sample ",Sample," for ArchR pre-processing step!"))
ArrowFiles <- paste0("ArrowFiles/",Sample,".arrow")

print("calculating doublet scores .....")
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "LSI", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  dimsToUse = 1:50,
  nTrials = 10,
  scaleDims = F,
  threads = 1,
  LSIParams = list(seed = 1,
  varFeatures = 50000,
  excludeChr = c("chrX", "chrY", "chrMT")),
  outDir = "002_doublets"
)
#Finally, an ArchR project can be generated.
print("creating ArchR project post doublet scores .....")

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = Sample[1], 	# generate an output directory before with this name
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)
##loading meta data
met_sams <- read.delim("Files/Ret_metdata.txt",row.names = 1)

proj$Stage <- met_sams[Sample[1],"Stage"]
proj$Age <- met_sams[Sample[1],"Age"]
proj$Sex <- met_sams[Sample[1],"Sex"]

##number of cells post intial QC
print("Number of cells before doublet filtering:")
nCells(proj)
plot.data <- data.frame(getCellColData(proj))
plot.data$barcode <- row.names(plot.data)
# Filtering putative doublets (pass 1).
# Here we want to be relatively lenient, so we select a ratio of 1 (i.e. filtering top 5% of barcodes for a sample of 5,000 cells).
# We are using Doublet enrichment as a cutoff.

proj <- filterDoublets(ArchRProj = proj,filterRatio = 1)

print("Number of cells post doublet filtering:")
nCells(proj)
##saving doublet cells
doublets <- plot.data$barcode[!(plot.data$barcode %in% getCellNames(proj))]
write(doublets, paste0("002_doublets/",Sample[1],"/",Sample[1],"_doublets.txt"))
##also figuring out high Fragment cells and cells with high Doublet enrichment
x <- mean(plot.data$nFrags)
z <- x+ 2*mad(plot.data$nFrags)
keep1 <- plot.data$nFrags > z
#remove doublet enrichment
keep2 <- plot.data$DoubletEnrichment > 4
put_doublets <- plot.data$barcode[keep1 | keep2]
length(put_doublets)
#total putative doublets that will be removed
cells_rem <- unique(c(put_doublets,doublets))
write(cells_rem,  paste0("002_doublets/",Sample[1],"/",Sample[1],"_rem_cells.txt"))
keep <- plot.data$barcode %in% cells_rem
plot.data <- plot.data[!keep,]
##subsetting object to filtered cells
proj <- proj[plot.data$barcode,]


#Now we can calculate our first iterative LSI. We are now extending our dimensions to 100 and the number of variable features to 100K.

proj <- addIterativeLSI(ArchRProj = proj,
                        useMatrix = "TileMatrix",
                        name = "IterativeLSI_qc",
                        iterations=5,
                        clusterParams = list(
                          resolution = c(0.1, 0.2, 0.4, 0.8),
                          sampleCells = 20000, ##this causes all the cells to be used as all samples are less than 29k
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



### post LSI ------

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

##extract LSI and GSM per sample
LSI <- getReducedDims(proj,"IterativeLSI")

save(LSI, file = paste0(Sample,"_LSI.RData"))

mdt <- getMatrixFromProject(proj, "GeneScoreMatrix")
genx <- data.frame(rowData(mdt))
mdt <- assay(mdt)
row.names(mdt) <- genx$name
save(mdt, file = paste0(Sample,"_GSM.RData"))
#saveArchRproj at last
print("Saving ArchR project:")
saveArchRProject(proj, outputDirectory = Sample[1])

q()
