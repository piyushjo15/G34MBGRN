## This is script for filtering peak matrix generated for Retina and generate
## robust peak set for downstream processing

suppressPackageStartupMessages({
  library(ArchR)
  library(tidyr)
  library(dplyr)
  library(Matrix)
})

##dirs
DIR_ATAC <- "Retina/ATAC/"
DIR_PEAK <- "Retina/ATAC/PeakANA/"
## obtaining a bed file of peak

########### filtering peaks###########

###checking if peaks are annotated ----------
load(paste0(DIR_PEAK,"PeakMatrixComRet.RData")) ##peak matrix
PeakMat <- assay(peak_mat)
region_names <- as.data.frame(peak_mat@rowRanges)
row.names(region_names) <- paste0(region_names$seqnames,":",region_names$start,"-",region_names$end)
head(region_names)
row.names(PeakMat) <- row.names(region_names)
# save(PeakMat,file = paste0(DIR_PEAK,"PeakMatrixComRetdata.RData") )

load(paste0(DIR_PEAK,"PeakSetComRet.RData"))##peakset
df <- data.frame(Peakset)
row.names(df) <- paste0(df$seqnames,":",df$start,"-",df$end)
head(df)
bed <- df[,c(1,2,3)]
head(bed)
bed$peaknames <- row.names(bed)
write.table(bed,paste0(DIR_PEAK,"combined_Retpeaks.bed"), sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

PS <- data.frame(Peakset) 
row.names(PS) <- paste0(PS$seqnames,":",PS$start,"-",PS$end)
head(PS)
##adding cluster info
cl <- PS$GroupReplicate
cl <-gsub("._."," ",cl,fixed = TRUE)
cl <- data.frame(X=cl)

cl2 <- cl %>% separate(X,c("A","B"), sep=" ")
head(cl2)
PS$Cluster <- cl2$A
head(PS)
rm(cl,cl2)
table(row.names(PS)==row.names(region_names))

###all peaks exist in peakset, and are in same order as in the PeakMatrix.!!!
## No filtering required for annotation!!

##Filtering for activity---------
load(paste0(DIR_ATAC,"plotdataRet_ATAC.RData"))
dim(plot.data)
clust <- sort(unique(PS$Cluster))
dim(PeakMat)
table(colnames(PeakMat) %in% row.names(plot.data)) ##all exists!
## rearranging the cell order for plotdata as thts easier
pd <- plot.data[colnames(PeakMat),]
table(colnames(PeakMat)==row.names(pd)) ##fixed!

##Ioannis's lapply function takes a lot of memory, so using for loop
peak_clust_f <- c()
for(id in clust){
  cell_i <- which(pd$predictedGroup==id)
  del  <- rowSums(PeakMat[,cell_i] > 0)/length(cell_i)
  peak_clust_f <- cbind(peak_clust_f,del)
  rm(del,cell_i)
}
rm(id)

colnames(peak_clust_f) <- clust
peak_clust_f[1:10,1:4]

save(PeakMat, PS, pd,peak_clust_f, file = paste0(DIR_PEAK,"peak_freq_cluster.RData"))
q()
##checking ----------------
load(paste0(DIR_PEAK,"peak_freq_cluster.RData"))

PS$max_freq <- apply(peak_clust_f, 1, max)
sum(PS$max_freq >= 0.01)
sum(PS$max_freq >= 0.03)
 sum(PS$max_freq >= 0.05)

ggplot(PS, aes(log10(max_freq))) +
  geom_histogram(bins = 30) +
  theme_classic()
p <- ggplot(PS, aes(x=log10(max_freq), y=log10(score))) +
  geom_hex() +
  theme_classic()
tiff(paste0(DIR_PEAK,"Peaks_maxFreq_byPeakScore.tiff"), unit = "in",  width=5, height=5, res = 300)
print(p)
dev.off()
p <- ggplot(PS, aes(x=as.factor(Reproducibility), y=log10(max_freq), fill=as.factor(Reproducibility))) +
  geom_violin() +
  geom_boxplot(notch = T, width=0.1) +
  scale_fill_brewer(palette = "Greens", name="") +
  geom_hline(yintercept = c(log10(0.03),log10(0.05)), color=c("red","orange"), lty="dashed") +
  xlab("Number of samples supporting peak") +
  theme_classic()

tiff(paste0(DIR_PEAK,"Peaks_maxFreq_byReproducibility.tiff"), unit = "in",  width=6, height=4, res = 300)
print(p)
dev.off()
## cut-off of 3% seems better as then most peaks are supported by 3 reps

cutoff <- 0.03

PS$robust <- PS$max_freq >= cutoff 

p <- ggplot(data = NULL, aes(x=seq(1, length(PS$max_freq), 20),y=sort(PS$max_freq[seq(1, length(PS$max_freq), 20)], decreasing = T), color=sort(PS$max_freq[seq(1, length(PS$max_freq), 20)], decreasing = T) >= cutoff )) +
  geom_point() +
  geom_hline(yintercept = cutoff , color="red", lty="dashed") +
  geom_vline(xintercept = which.min(sort(PS$max_freq, decreasing = T) >= cutoff ), color="red", lty="dashed") +
  scale_color_manual(values = c("gray70", "deepskyblue3"), guide=F) +
  annotate(geom="text",x=which.min(sort(PS$max_freq, decreasing = T) >= cutoff ) - 5e4, y=0.2, label=paste0("Robust peaks\nFreq >= ",cutoff ,"\nN=", which.min(sort(PS$max_freq, decreasing = T) >= cutoff )), color="deepskyblue3") +
  ylab("Maximum fraction of active cells in a cluster") +
  xlab("Peaks (decreasing maximum activity)") +
  theme_classic()


tiff(paste0(DIR_PEAK,"Peaks_robustPeak_identification.tiff"), unit = "in",  width=5, height=5, res = 300)
print(p)
dev.off()
save(PeakMat, PS, pd, file = paste0(DIR_PEAK,"Peaks_robust.RData"))
##write robust and normal peaks into bed file
bed <- PS[,c("seqnames",  "start",    "end")]
head(bed)
write.table(bed, file=paste0(DIR_PEAK,"combined_Retpeaks.bed"),
            row.names = FALSE, col.names = FALSE,quote = FALSE, sep="\t")
bed <- PS[PS$robust,c("seqnames",  "start",    "end")]
head(bed)
write.table(bed, file=paste0(DIR_PEAK,"Retina_allPeaks_robust.bed"),
            row.names = FALSE,col.names = FALSE, quote = FALSE, sep="\t")
##these peaks are sorted later

q()
##this generated Retina_allPeaks_robust.bed was sorted and filtered to remove blacklisted regions
#sort -k 1,1 -k2,2n Retina_allPeaks_robust.bed > Retina_allPeaks_robust.sorted.bed
#bedtools subtract -A -a Retina_allPeaks_robust.sorted.bed -b Files/hg38-blacklist.v2_noncanchr.bed > Retina_allPeaks_robust.fil.sorted.bed
