## This script performs peak calling on integrated data

############# Peak Calling Integrated Data ##########################
#PEAK CALLING

suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(hexbin)
  library(BSgenome.Hsapiens.UCSC.hg38) 
})



set.seed(456)
addArchRThreads(threads = 12)
#setwd
#Next, I load the reference genome
addArchRGenome("hg38")
Sample <- "CombinedATAC_Ret"
DIR="Retina/ATAC/"
# load the data
proj <- loadArchRProject(path = Sample)


#Creating PseudoBulk Replicates
proj <- addGroupCoverages(
  ArchRProj = proj,
  sampleRatio = 0.8, 
  returnGroups = F,
  force = T,
  groupBy = "predicted" ,
  minCells = 100,
  maxCells = 1000,
  minReplicates = 2,
  maxReplicates = 8,
  maxFragments = 50 * 10^6,
  useLabels= TRUE
)

saveArchRProject(ArchRProj = proj, load = FALSE)

###### 4. Peak Calling ##########

#Define path to MACS2
pathToMacs2 <- findMacs2()

#Peak Calling
proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "predicted",
  pathToMacs2 = pathToMacs2,
  reproducibility = "2",
  excludeChr = c("chrY", "chrMT"), 
  method = "q",
  cutOff = 0.01, 
  extendSummits = 250,
  force = T
)
print(paste0("Finished PeakCalling for Sample:", Sample))


#Add PeakMatrix
proj <- addPeakMatrix(proj,
                      ceiling = 5,
                      binarize = F,
                      force = T)

####### Save ArchR project ##########
saveArchRProject(ArchRProj = proj, load = FALSE)

##obtain peak matrix
Peakset <- getPeakSet(proj)
save(Peakset, file = paste0(DIR,"PeakANA/PeakSetComRet.RData"))
peak_mat <- getMatrixFromProject(proj, "PeakMatrix")
save(peak_mat, file = paste0(DIR,"PeakANA/PeakMatrixComRet.RData"))

q()
