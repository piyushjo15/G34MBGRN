## This script is used for generating  combined ATAC ArchR project
## This is for copying original unaltered arrows files into a new folder
## for integrated analysis of all samples together.
suppressPackageStartupMessages({
  library(ArchR)
})


addArchRThreads(threads =6)
addArchRGenome("hg38")
Sample <- "CombinedATAC"

DIR="ATAC/"
#Define Arrow files
Arrow_Path <- paste0(DIR,"ArroWFiles/")
ArrowFiles <- list.files(Arrow_Path)
ArrowFiles <- paste0(Arrow_Path,ArrowFiles)

proj <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = Sample,
    copyArrows = TRUE
)

#saveArchRproj in between
print("Saving ArchR project:")
saveArchRProject(proj, outputDirectory = Sample)

## removing cells that did not pass QC filtering in each sample 
plot.data <- load("ATAC/fixdplotdata.RData")
keep <- plot.data$ATAC_cells

proj_com <- subsetArchRProject(proj_com,cells = keep,outputDirectory = Sample, force = TRUE)
print("Filtered cells")

proj_com <- loadArchRProject(Sample)
### obtaining joint representation of all the ATAC data

######################## LSI ################################
#calculate iterative LSI
proj_com <- addIterativeLSI(ArchRProj = proj_com,
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
proj_com <- addClusters(input = proj_com,
                    name = "Clusters_int",
                    reducedDims = "IterativeLSI_int",
                    method = "Seurat",
                    force = T,
                    resolution=1,
                    corCutOff = 0.75,
                    scaleDims = FALSE,
                    seed = 1)


# UMAp calculation
proj_com <- addUMAP(ArchRProj = proj_com,
                name = "UMAP_int",
                reducedDims = "IterativeLSI_int",
                minDist = 0.2,
                metric = "cosine",
                nNeighbors = 25,
                force = T,
                seed = 1,
                scaleDims = F,
                corCutOff = 0.75)
#saveArchRproj
print("Saving ArchR project:")
saveArchRProject(proj_com, outputDirectory = Sample)

print("Finished Combining Arrow Files to one proj")

## Adding gene expression data to each cell in this combined ATAC ArchR project
## This will be used for gene to peak links analysis
MB_Samples <- unique(proj_com$Sample)

GEX_list <- list()
for(i in 1:length(MB_Samples)){
  RNA_DIR <- paste0(DIR,MB_Samples[i],"/Imputation/")
  load(paste0(RNA_DIR,MB_Samples[i],"_imputed.RData"))
  x <- mdt_x
  GEX_list[[MB_Samples[i]]] <- x
}

#cbind all matrices to have 1 matrix!
combined_matrix <- do.call(cbind, GEX_list)
dim(combined_matrix) 

#subset for cells present in proj
cells <- getCellNames(proj_com)
combined_matrix <- combined_matrix[, cells]
dim(combined_matrix)

save(combined_matrix,file="Multiome/combined_GEX_matrix.RData")

#save the combined matrix as SE
gr_obj <- rtracklayer::import("/b06x-isi/b062/m-o/MBATAC/data/t716r/RNA/gencdH38p13r37CR_genesN.filtered.bed")
SE <- SummarizedExperiment(assays = list(logcounts=combined_matrix), rowRanges = gr_obj)

#add GeneExpressionMatrix
proj_com <- addGeneExpressionMatrix(
  input = proj_com,
  seRNA = SE,
  chromSizes = getChromSizes(proj_com),
  excludeChr = c("chrM", "chrY")
)

saveArchRProject(proj_com)
save(SE,file="Multiome/combined_GEX_matrix_SE.RData")



