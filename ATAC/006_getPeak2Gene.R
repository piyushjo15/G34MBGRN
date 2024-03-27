## This script calculates peak to gene links in combined ATAc ArchR object
## post peak calling based on clusters identified from GRN analysis

#P2G links ArchR

suppressPackageStartupMessages({
  library(ArchR)
})

set.seed(456)
addArchRThreads(threads = 12)

addArchRGenome("hg38")
Sample <- "CombinedATAC_fil"
DIR="ATAC/"
setwd( paste0(DIR,Sample,"/"))

proj_com <- loadArchRProject(path = Sample)


#add Peak to gene links
proj_com<- addPeak2GeneLinks(
  ArchRProj = proj_com,
  reducedDims = "IterativeLSI_int",
  useMatrix = "GeneExpressionMatrix"
)


p2g <- getPeak2GeneLinks(
  ArchRProj = proj_com,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE #FALSE, df is returned
)

p2g_granges <- getPeak2GeneLinks(
  ArchRProj = proj_com,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = TRUE #TRUE, granges returned
)

write.table(p2g,"Peak2Gene_df.txt",col.names = TRUE,sep="\t")
write.table(p2g_granges,"Peak2Gene_Granges.bed",col.names = T,sep = "\t")


p2geneDF <- metadata(proj_com@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), ":", start(.), "_", end(.))})[p2geneDF$idxATAC]

write.table(p2geneDF,paste0(Sample,"_Peak2Gene_df_withNames.txt"),sep = "\t",col.names = TRUE,row.names = F)

#save
saveArchRProject(proj_com)






