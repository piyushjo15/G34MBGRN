## This script adds missing links of TF-region-target that are missed by SCENIC+ 
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
#get a TF to region matrix
args <- commandArgs(TRUE)
Sample <- args[1]

DIR = "~Retina/SCENIC/scenicplus/"
#read metadata
setwd(paste0(DIR, Sample, "/output/"))

df <- read.csv(paste0(Sample,"_eRegulon_metadata.csv"), sep = ";")
df2 <- df[df$is_extended=="True",]
gg <- as.character(df2$Gene_signature_name)
ggx <- grep("+_+",gg, fixed = TRUE)

##extracting enhancer identified in this sample
df3 <- df2[ggx,]
tfs <- unique(df3$TF)

##unique enhancer regions
reg_adj <- data.frame(Region=df3$Region)
reg_adj <- reg_adj %>% separate(Region,c("A","B","C"))

keep <- duplicated(df3$Region)
table(keep)
reg_adj <- reg_adj[!keep,]

head(reg_adj)
##this bed needs to be sorted
write.table(reg_adj, file=paste0(Sample,"_enhancers.bed"), 
            sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)

################# Processing to fix missing links ####################
###extracting regions associated with a TF for positive regulations
TF_region <- list()
for (x in tfs){
  del <- df3[df3$TF==x,c("TF","Region")]
  del <- del[!duplicated(del$Region),] ##remove duplicates
  TF_region[[x]] <- del
  rm(del)
}
rm(x)
#get a region to gene matrix
regs <- unique(df3$Region)
region_tar <- list()
for (x in regs){
  del <- df3[df3$Region==x,c("Region","Gene", "R2G_importance", "R2G_rho", "R2G_importance_x_rho")]
  XX <- paste0(del$Region,"_",del$Gene)
  del <- del[!duplicated(XX),]#remove duplicates
  region_tar[[x]] <- del
  rm(del,XX)
}
rm(x)

##now I want three things, TF, regions and gene links
GRNS <- c()
for(x in tfs){
  regions <- as.character(TF_region[[x]]$Region)
  targets <- c()
  for (y in regions){
    del <- region_tar[[y]]
    targets <- rbind(targets,del)
    rm(del)
  }
  targets$TF <- x
  GRNS <- rbind(GRNS,targets)
  rm(y,targets,regions)
}
rm(x)
table(duplicated(paste0(GRNS$Region,"_",GRNS$Gene,"_",GRNS$TF))) ##checking duplicated connections, ALL false

row.names(GRNS) <- seq(dim(GRNS)[1])
head(GRNS)
#get a region to gene matrix
fp <- "~/MBsnANA/Retina/SCENIC/" 
TF2G <- read.delim(paste0(fp,Sample,"/cor",Sample,"_cleaned.tsv")) ### fix the file loc

colnames(TF2G)[c(3,5)] <- c("TF2G_importance","TF2G_rho")
##selecting only pos corr
TF2G <- TF2G[TF2G$TF2G_rho>=0.05,] ##adding cut-off rho
GRNS2 <- c()
for(x in tfs){
  TF2G_sel <- TF2G[TF2G$TF==x,]
  ##add TF-TF linksif not present
  if(!(x %in% TF2G_sel$target )){
    imp <- max(TF2G_sel$TF2G_importance)+ 1e-05
    TF2G_sel <- rbind(TF2G_sel, c(x,x,imp,1,1))
  }
  tars1 <- unique(TF2G_sel$target)
  GRNS_sel <- GRNS[GRNS$TF==x,]
  tars2 <- unique(GRNS_sel$Gene)
  targets <- intersect(tars1,tars2)
  ##subset to shared targets
  GRNS_sel <- GRNS_sel[GRNS_sel$Gene %in% targets,]
  GRNS_sel$TF2G_importance <- GRNS_sel$TF2G_rho <- NA
  
  for (y in targets){
    GRNS_sel[GRNS_sel$Gene==y,"TF2G_importance"] <-TF2G_sel[TF2G_sel$target==y,"TF2G_importance"] 
    GRNS_sel[GRNS_sel$Gene==y,"TF2G_rho"] <-TF2G_sel[TF2G_sel$target==y,"TF2G_rho"]
  }
  GRNS2 <- rbind(GRNS2,GRNS_sel)
  rm(y,targets,GRNS_sel,TF2G_sel,tars1, tars2 )
}
rm(x)
GRNS2 <- GRNS2[,c(6,1:5,7,8)]
row.names(GRNS2) <- seq(dim(GRNS2)[1])
tail(GRNS2)
write.table(GRNS2,paste0(Sample,"_eRegulons_fix.txt"),quote = FALSE, row.names = FALSE, sep ="\t")

##getting TF-GRN
def <- GRNS2[,c("TF","Gene")]
tfs <- unique(def$TF)
regs <- list()
for(x in tfs){
  defx <-def[def$TF==x,]
  regs[[x]] <- unique(c(defx$Gene,x))
  rm(defx)
}
rm(x)
save(regs, file = paste0(Sample,"_gmt.RData"))
#make a bed file
GRNS2 <- GRNS2 %>% separate(Region,c("chr","start","end"))
head(GRNS2)
GRNS2$ID <- paste0(GRNS2$TF,"_",GRNS2$start,"_",GRNS2$Gene)
bed <- GRNS2[,c("chr","start","end","ID")]
head(bed)
write.table(bed,paste0(Sample,"_enhancersv2.bed"),quote = FALSE, row.names = FALSE, col.names = FALSE, sep ="\t")

##extracting links
REGS <- GRNS2[,c("TF","chr","start","end","ID","Gene")]
save(REGS,file = paste0(Sample,"_TF2R2G.RData"))

q()