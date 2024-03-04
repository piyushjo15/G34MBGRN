## This script obtained regions associated with Tf-target relation for important GRNs
## across samples
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
#Dirs
DIR_BED <- "Retina/ATAC/SCENIC/BEDS/"
DIR_GRN <- "Retina/ATAC/SCENIC/scenicplus/"

##
Sample = commandArgs(trailingOnly=TRUE)

##get the GRNs
load(paste0(DIR_GRN,"TFselcl_GRN_Retina.RData")) #G34MB GRNs, loading gene sets identifed from combine GRNS
tfs <- names(regs) ##these are the candidate TFs
#write(tfs, file =paste0(DIR_BED,"final_tfs.txt"))
##load sample associated GRNs
GRN <- read.delim(paste0(DIR_GRN,Sample[1],"/output/",Sample[1],"_eRegulons_fix.txt"))

##for each sample get all the regions associated with TF-target from combined GRNs
keep <- GRN$TF %in% tfs
table(keep)
GRN <- GRN[keep,] ##subset to interesting TF
head(GRN)
## now collecting for each TF its associated target as present in consistent GRN 
## and region associated wih it
tfs2 <- unique(GRN$TF)
head(tfs2)
for(x in tfs2){
  GRNx <- GRN[GRN$TF==x,]
  tars <- regs[[x]] ##get associated targets
  GRNy <- GRNx[GRNx$Gene %in% tars,] ##subset
  if(dim(GRNy)[1]>=1){
    GRNz <- GRNy %>% separate(Region,c("chr","start","end"))
    GRNz$Link <- paste0(GRNz$TF,"_",GRNz$start,"_",GRNz$Gene)
    bed <- GRNz[,c("chr","start","end","Link")]
    #print("wringing bed ...")
    write.table(bed, file = paste0(DIR_BED,x,"_",Sample[1],"_GRN.bed"),quote = FALSE, row.names = FALSE, col.names = FALSE, sep= "\t")
    rm(GRNy, GRNz, bed)
  }
  rm(GRNx, tars)
}


q()