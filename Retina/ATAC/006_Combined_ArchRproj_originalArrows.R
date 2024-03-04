# This script is used for a combined project
## This is for copying original unaltered arrows files into a new folder
## for integrated analysis of all samples together.
suppressPackageStartupMessages({
  library(ArchR)
  library(hexbin)
  library(BSgenome.Hsapiens.UCSC.hg38) 
})


addArchRThreads(threads =1)
#Next, I load the reference genome
addArchRGenome("hg38")
Sample <- "CombinedATAC_Ret"

DIR="Retina/ATAC/"

load(paste0(DIR, "com_RNA_ATAC_pd_scJ.RData"))
pd_ATAC <- plot.data[plot.data$Platform=="ATAC",]

#Define Arrow files
Arrow_Path <- paste0(DIR,"original_processed_AF/")
ArrowFiles <- list.files(Arrow_Path)
ArrowFiles <- paste0(Arrow_Path,ArrowFiles)


proj <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = Sample,
    copyArrows = TRUE
)
cells <- getCellColData(proj)
table(row.names(cells) %in% row.names(pd_ATAC))
proj <- proj[,row.names(pd_ATAC)]
proj@cellColData <-DataFrame(pd_ATAC)


######################## LSI ################################
#calculate iterative LSI
proj <- addIterativeLSI(ArchRProj = proj,
                        useMatrix = "TileMatrix",
                        name = "IterativeLSI_int",
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
                    name = "Clusters_int",
                    reducedDims = "IterativeLSI_int",
                    method = "Seurat",
                    force = T,
                    resolution=1,
                    corCutOff = 0.75,
                    scaleDims = FALSE,
                    seed = 1)


# UMAP
proj <- addUMAP(ArchRProj = proj,
                name = "UMAP_int",
                reducedDims = "IterativeLSI_int",
                minDist = 0.2,
                metric = "cosine",
                nNeighbors = 25,
                force = T,
                seed = 1,
                scaleDims = F,
                corCutOff = 0.75)
#saveArchRproj in between
print("Saving ArchR project:")
saveArchRProject(proj, outputDirectory = Sample)




q()