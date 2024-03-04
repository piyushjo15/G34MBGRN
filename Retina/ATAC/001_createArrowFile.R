

#This script is used to generate Arrow files from fragment files

suppressPackageStartupMessages({
  library(ArchR)
  library(tibble)
  library(hexbin)
  library(BSgenome.Hsapiens.UCSC.hg38) 
})

set.seed(456)
addArchRThreads(threads =12) 
#define the Sample ID
Sample = commandSample(trailingOnly=TRUE)


#Next, I load the reference genome
addArchRGenome("hg38")
load("Files/GENCDH38p13r37_ann_4_archr.RData")

#define the input file
inputFiles <- paste0("Retina/",Sample,"/atac_fragments.tsv.gz")

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
