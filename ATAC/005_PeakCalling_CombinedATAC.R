## This script calls peak on combined ATAC ArchR object

############# Peak Calling Integrated Data ##########################
#PEAK CALLING

suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38) 
})



set.seed(456)
addArchRThreads(threads = 12)
addArchRGenome("hg38")
Sample <- "CombinedATAC"

DIR_ATAC <- "ATAC/"
DIR_GRN <- "SCENIC/GRNana/"

setwd( paste0(DIR,Sample,"/"))

# load the data
proj_com <- loadArchRProject(path = Sample)

plot.data_withNormal <- getCellColData(proj_com)

#### load the plot.data of integrated data ###############################
##This metadata is post GRN analysis. Clusters are identified from GRN analysis.
## Only tumor cells are included
plot.data <- load("SCENIC/GRNana/plotdata_G34MBaucKNN.RData") 

#subset for cells where ATAC is present
plot.data_ATAC <- plot.data[plot.data$ATAC == "Yes",]
nCells(proj_com)
nrow(plot.data_ATAC) 

#subset for cells present in plot.data (= normal cells are removed)

keep <- plot.data_ATAC$ATAC_cells
cells <- getCellNames(proj_com)
#check if all keep cells are present in cells!
keep <- cells[cells %in% keep] 

#subset for these cells (normal cells are removed!)
subsetArchRProject(proj_com,
                   cells = keep,
                   outputDirectory = paste0(Sample,"_fil/"),
                   force=TRUE)

proj_com <- loadArchRProject(Sample, outputDirectory=paste0(Sample,"_fil/"))

## Adding tumor cell metadata to this filtered combined ATAC ArchR object
colData <- getCellColData(proj_com)
cellsATAC <- rownames(colData)
cells_plotdata <- plot.data_ATAC$ATAC_cells

# Get the order 
order <- match(cellsATAC, cells_plotdata)
# Reorder vector 1
plot.data_ATAC <- plot.data_ATAC[order,]

#since the order of the barcode is the same I can directly add the columns prom plot.data to the ArchR proj

add <- plot.data_ATAC[,c("Axis","Annotation","Batch","ATAC_Batch","Group","Subtype","Age","Sex","Imputed","KNN_cl","nfUMAP1","nfUMAP2")]
colData <- cbind(colData,add)
proj_com@cellColData <- colData

########## add the UMAP ########################################

#add UMAP
UMAP <- data.frame(plot.data_ATAC$nfUMAP1, plot.data_ATAC$nfUMAP2)
colnames(UMAP) <- c("LSI#UMAP_Dimension_1","LSI#UMAP_Dimension_1")
rownames(UMAP) <- plot.data_ATAC$ATAC_cells
proj_com@embeddings@listData$UMAP$df <- UMAP

#Creating PseudoBulk Replicates
proj_com <- addGroupCoverages(
  ArchRProj = proj_com,
  sampleRatio = 0.8, 
  returnGroups = F,
  force = T,
  groupBy = "KNN_cl" , ## this 
  minCells = 100,
  maxCells = 1000,
  minReplicates = 5,
  maxReplicates = 15,
  maxFragments = 50 * 10^6,
  useLabels= TRUE
)


###### Peak Calling ##########

#Define path to MACS2
pathToMacs2 <- findMacs2()

#Peak Calling
proj_com <- addReproduciblePeakSet(
  ArchRProj = proj_com, 
  groupBy = "KNN_cl",
  pathToMacs2 = pathToMacs2,
  reproducibility = "2",
  excludeChr = c("chrY", "chrMT"), 
  method = "q",
  cutOff = 0.01, 
  extendSummits = 250,
  force = T
)

#Adding Peak matrix
proj_com <- addPeakMatrix(proj_com,
                      ceiling = 5,
                      binarize = F,
                      force = T)


##########ChromVAR Deviations Enrichment############
library(chromVAR)
library(chromVARmotifs)
#add a set of background peaks used in computing deviations
proj_com <- addBgdPeaks(proj_com, force = T)

#add Jaspar annotation
load("Files/JASPAR2022PWMlist.RData") #pwmList and TF
names(pwmList) <- TF$TF 

proj_com <- addMotifAnnotations(ArchRProj = proj, motifSet = "JASPAR2020", annoName = "jaspar", force = T, motifPWMs = pwmList)

#compute per-cell deviations across all motif annotations -> Matrix is generated

proj_com <- addDeviationsMatrix(
  ArchRProj = proj_com, 
  peakAnnotation = "jaspar",
  matrixName = "MotifMatrix", #was added as MotiMatrix!!!! the mistake was changed after matrix was added!
  force = TRUE,
  threads=1
)

saveArchRProject(ArchRProj = proj_com, load = FALSE)

print(paste0("Finished PeakCalling for Sample:", Sample))
