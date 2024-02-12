
##samples used for SCENIC+ analysis
sams <- read.delim("Files/MBsamples.txt", header = FALSE)
sams <- sams$V1
DIR_GRNana <- "SCENIC/GRNana/"
DIR_GRN <- "SCENIC/scenicplus/"


## Using obtained marker GRNs per sample, identify important TF-GRNS
com <- c()
for(x in sams){
  load(paste0(DIR_GRNana,x,"_top3GRN_AUC.RData")) 
  rm(ups)
  com <- c(com,all_reg)
  rm(all_reg)
}
rm(x)
df <- data.frame(table(com))
df <- df[df$Freq>1,] ## remove sample specific GRNs, GRN is a cluster marker in atleast 2 sample
regs_com <- as.character(df$com)

#regs_com <- sort(unique(com))
## The GRN can exist in a sample without being cluster specific for example MYC or CRX in Group 3
## Since that GRN is still playing import role in that specific sample, we will use all samples where
## a GRN was identified to get a list of consistent target genes
ADJ <- list()
## gathering imp GRNs from each sample where it was identified
for(x in sams){
  z <- gsub("MSN","MSM",x)
  if(x=="SN407"){z="SA194"}
  load(paste0(DIR_GRN,z,"/output/metadata/",x,"_gmt.RData")) ##get identified GRNs component per sample
  regx <- names(regs)
  regx <- intersect(regx,regs_com) ##intersect for only imp GRns
  GS <- list()
  ## Get a list of TF-targets for a TF from each sample
  for(y in regx){
    GS[[y]] <- regs[[y]]
  }
  ADJ[[x]] <- GS ## subsetted list of important GRNs per sample
  
  rm(GS, regx,y, regs,z)
}
rm(x)


## Now I want to get a list of targets for each TF that appeared in different samples,
## and create a GRN that is sort of consistent
regs <- list()
for( y in regs_com){
  sel <- c()
  for(x in sams){
    del <- ADJ[[x]]
    tfc <- names(del)
    if(y %in% tfc) {
      del2 <- del[[y]]
      sel <- c(sel,del2)
      rm(del2)
    }
    rm(del, tfc)
  }
  rm(x)
  
  if(length(sel)>1){
    ##removing lowly represented genes
    tb <- data.frame(table(sel))
    # taking only those targets that appeared at least 20% of time (samples) 
    # the GRN was detected, or atleast thrice
    co <- max(round(max(tb$Freq)/5)-1,2) 
    keep3 <- tb$Freq>co
    tb <- tb[keep3,]
    regs[[y]] <- as.character(tb$sel)
    rm(keep3)
  }
  rm(sel,tb)
}
rm(y)

for(y in regs_com){
  del <- length(regs[[y]])
  if(del<15){regs[[y]] <- NULL} # remove gene-sets having lower than 15 genes
}
rm(y)
lengths(regs)
save(regs, file = paste0(DIR_GRNana,"TFselcl_GRN_AUC_G34MB.RData"))
##finding jaccard and overlap index for identified GRN sets
library(RColorBrewer)
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
overlap <- function(a, b) {
  intersection = length(intersect(a, b))
  min_size = min(length(a),length(b))
  return (intersection/min_size)
}

sams <- names(regs)
JS <- matrix(0, nrow = length(sams),ncol = length(sams))
OS <- matrix(0, nrow = length(sams),ncol = length(sams))

for(i in 1:length(sams)){
  for (j in 1:length(sams)){
    JS[i,j] <- jaccard(regs[[i]],regs[[j]])
    OS[i,j] <- overlap(regs[[i]],regs[[j]])
  }
}
colnames(JS) <- row.names(JS) <-
  colnames(OS) <- row.names(OS) <- sams
sel <- read.delim("GRNs.txt")
sel <- sel$GRNs
mc <- colorRampPalette(c("white",brewer.pal(9,"BuPu")))(100)
mb <- ((seq(1:101)-1)/100)
hp <- pheatmap::pheatmap(JS[sel,sel], 
                         cluster_rows = FALSE, 
                         cluster_col = FALSE,
                         fontsize = 6,show_colnames = FALSE,
                         border_color = NA,
                         breaks = mb, color = mc)
tiff(paste0(DIR_GRNana,"TFGRN_JS.tiff"), width = 8, height = 8.5, units = "in", res = 300)
print(hp)
dev.off()
mc2 <- colorRampPalette(c("white",brewer.pal(9,"Greys")))(100)
hp <- pheatmap::pheatmap(OS[sel,sel], 
                         cluster_rows = FALSE, 
                         cluster_col = FALSE,
                         fontsize = 6,show_colnames = FALSE,
                         border_color = NA,
                         breaks = mb, color = mc2)
tiff(paste0(DIR_GRNana,"TFGRN_OS.tiff"), width = 8, height = 8.5, units = "in", res = 300)
print(hp)
dev.off()
q()