suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(data.table)
  library(Matrix)
  library(SummarizedExperiment)
})
##dirs
DIR_ATAC <- "Retina/ATAC/"
DIR_PEAK <- "Retina/ATAC/PeakANA/"
DIR_GRN <- "Retina/ATAC/SCENIC/scenicplus/"

#
###analysis ------
##loading peak matrix
load(paste0(DIR_PEAK,"Peaks_robust.RData"))
load(paste0(DIR_GRN,"pdRNA_ATAC_fMNN_AUC.RData"))
pd <- plot.data[plot.data$Platform=="ATAC",c("predicted","Sample")]
rm(plot.data, rd , rd_cor)
bed <- read.delim(paste0(DIR_PEAK,"Retina_allPeaks_robust.fil.sorted.bed"), header = FALSE)
sel <- paste0(bed$V1,":",bed$V2,"-",bed$V3)
head(sel)
##pd is plotdata, already rearranged by cellids as present in PeakMatrix.M
##filter to only robust peaks and removing blacklisted region
PS$Final <- "No"
keep <- row.names(PS) %in% sel
table(keep)
PS[keep,"Final"] <- "Yes"

keep <- row.names(PeakMat) %in% sel
table(keep)
PeakMat_sel <- PeakMat[keep,row.names(pd)]
rm(PeakMat)
dim(PeakMat_sel)
table(row.names(pd)==colnames(PeakMat_sel))
##grouping by cluster
cell_counts <- group_by(pd, predicted) %>%
  ## Keeping only major lineages
  dplyr::count() %>%
  mutate(sample=(predicted))
cell_counts.filtered <- filter(cell_counts, n >=50) 
print(cell_counts.filtered$sample)

##this is by predicted annotation
pseudobulks <- do.call(cbind, parallel::mclapply(1:nrow(cell_counts.filtered), function(i){
  cells <- row.names(pd)[pd$predicted==cell_counts.filtered$predicted[i]]
  pseudo <- rowSums(PeakMat_sel[, cells])
  return(pseudo)
}, mc.cores = 10))

dim(pseudobulks)
colnames(pseudobulks) <- unique(cell_counts.filtered$sample)
pseudobulks[1:10,1:10]
pseudobulks2 <- t(t(pseudobulks)/colSums(pseudobulks) * 1e4) ##not doing 1e6 like Ioannis

##Ioannis's normalization
max_cpm <- apply(pseudobulks2, 1, max)
#hist(log10(max_cpm))

#Standardising: fraction of max cpm value

pseudobulks.std <- pseudobulks2/max_cpm

pseudobulks.std[1:5, 1:5]
###logtrasnformation, not following Ioannis's CPM normalization
mat <- log2(pseudobulks2 + 1)
save(pseudobulks,pseudobulks.std, pd,mat, file = paste0(DIR_PEAK,"Pseudoblk_peaks_bypredicted.RData"))
q()