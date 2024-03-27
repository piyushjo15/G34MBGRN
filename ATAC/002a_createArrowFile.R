
## This script is used to generate Arrow files from fragment files obtained from 
## cellranger-atac/arc alignment

suppressPackageStartupMessages({
  library(ArchR)
  library(tibble)
  library(hexbin)
  library(BSgenome.Hsapiens.UCSC.hg38) 
})

set.seed(456)
addArchRThreads(threads =12) 
Sample = commandSample(trailingOnly=TRUE)

#Next, I load the reference genome
addArchRGenome("hg38")
DIR="ATAC/ArrowFiles/"

#custom genome annotation based on gencode genome version used for RNA alignment
load("Files/GENCDH38p13r37_ann_4_archr.RData") 

#location of fragment file from cellranger alignment
inputFiles <- paste0("ATAC_alignment/",Sample,"/outs/atac_fragments.tsv.gz")

#generate ArrowFiles
setwd(DIR)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = Sample,
  geneAnnotation=new_gen_ann,
  minTSS = 3,  
  minFrags = 3000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  QCDir = "001_barcode_qc",
  promoterRegion = c(2000,100),
  force = T
)
warning()
print("Generated Arrow File ")


q()
