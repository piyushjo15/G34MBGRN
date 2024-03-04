## Using peaks overlapping reglon of identified regulons and using PeakMatrix
## to identify count activity
suppressPackageStartupMessages({
  library(Matrix)
  library(pheatmap)
  library(RColorBrewer)
})

##dirs
DIR_ATAC <- "Retina/ATAC/ArchR/"
DIR_PEAK <- "Retina/ATAC/PeakANA/"
DIR_GRN <- "Retina/ATAC/SCENIC/scenicplus/"
DIR_GRNx <- "SCENIC/GRNana/"

#####
##loading peak matrix
load(paste0(DIR_PEAK,"Pseudoblk_peaks_bypredicted.RData")) ##already Pseudobulks by predicted cell-type
##TFs
tfs <- readLines("Files/GRNs.txt")
##for each TF I will identify the associated peaks and find total reads associated with them
PSB_REGS <- c()
for(x in tfs){
  peaks <- read.delim(paste0(DIR_PEAK,x,"_final_peaks.bed"), header = FALSE)
  sel <- paste0(peaks$V1,":",peaks$V2,"-",peaks$V3)
  keep <- (sel %in% row.names(pseudobulks))
  sel <- sel[keep]
  psb_sel <- colSums(pseudobulks[sel,])
  PSB_REGS <- rbind(PSB_REGS,psb_sel)
  rm(peaks,sel,keep,psb_sel)
}
rm(x)
row.names(PSB_REGS) <- tfs
PSB_REGS[1:4,1:4]

pseudobulks2 <- t(t(PSB_REGS)/colSums(PSB_REGS) * 1e4) ##not doing 1e6 like Ioannis

mat <- log2(pseudobulks2 + 1)
##for ordering cluster
cl_lv <- c("Pro_early","Pro_late1","Pro_late2","Pro_late_II1","Pro_late_II2",
           "RGC_e1","RGC_e2","RGC","AC/HC_pro","AC_e","AC_1","AC_2",
           "HC_e1","HC_e2","HC",
           "Cones_e1","Cones_e2","Cones_l",
           "Rod_e","Rods_l","Rods_II",
           "BiP_e1","BiP_e2","BiP_on1","BiP_on2","BiP_on3","BiP_off1","BiP_off2","BiP_off3",
           "Astro","MG")
cl_lv <- cl_lv[cl_lv %in% colnames(pseudobulks)]
zscore_c <- function(x){
  sdx <- sd(x)
  mx <- mean(x)
  ot <- boxplot.stats(x)$out
  q <- boxplot.stats(x)$stats
  keep <- (x < q[1])
  x[keep] <- q[1]
  keep <- (x > q[5])
  x[keep] <- q[5]
  x2 <- x-mx
  z <- x2/sdx
  return(z)
}
mat2 <- t(apply(mat, 1, zscore_c))
mat2[mat2>2] <- 2
mat2[mat2<(-2)] <- (-2)
colnames(mat2) <- colnames(mat)
##heatmap
mc <- colorRampPalette(c(rev(brewer.pal(9,"Purples")),"white",brewer.pal(9,"OrRd")))(100)
mc <- colorRampPalette(c("white",brewer.pal(9,"Purples")))(100)
mc2 <- colorRampPalette(c("steelblue3","grey98","orangered3"))(100)


hp <- pheatmap(mat2[tfs,cl_lv],
        cluster_rows =  F,
        cluster_cols = F,
        color = mc2,
        border = NA)
tiff(paste0(DIR_GRN,"G34MB_enhancer_fixed_Retina.tiff"), width = 12, height = 18, units = "in", res = 300)
print(hp)
dev.off()
