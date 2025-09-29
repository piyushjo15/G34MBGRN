## This script annotates summits associated with EOMES and CRX ChIP
## and find overlapping target genes
suppressPackageStartupMessages({
  library(rGREAT)
  library(ggplot2)
  library(gridExtra)
  library(GenomicRanges)
  library(tidyverse)
})

DIR_OUT <- "/path/to/output/"
setwd(DIR_OUT)
## 1. GREAT annotation ------
## For peaks assocaited with each TF
eomes <- read.delim("macs3/EOMES.com.bed",header = FALSE)
crx <- read.delim("macs3/CRX.com.bed",header = FALSE)
eomes$ID <- "EOMES"
crx$ID <- "CRX"
beds <- rbind(eomes,crx)
head(beds)
beds$V4 <- paste0(beds$V1,":",beds$V2,"-",beds$V3)
beds$strand <- "."

## duplicated peaks
table(duplicated(beds$V4))
cre <- beds$V4[duplicated(beds$V4)]
keep <- beds$V4 %in% cre
beds2 <- beds[keep,]
beds2$ID <- "Both"

beds2 <- beds2[duplicated(beds2$V4),]
beds <- beds[!keep,]
beds <- rbind(beds,beds2)
row.names(beds) <-beds$V4
rm(beds2)
table(beds$ID)

## covnert to GRanges object
PS_gr <- makeGRangesFromDataFrame(beds,
                                  keep.extra.columns=FALSE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field=c("V1"),
                                  start.field="V2",
                                  end.field="V3",
                                  strand.field="strand",
                                  starts.in.df.are.0based=TRUE)

##dummy geneset to avoid GO analysis
gs <- list(What=c("ENSG00000186092"))
res <- great(PS_gr,gs,"Gencode_v37",min_gene_set_size = 1)

R2G <- getRegionGeneAssociations(res)
R2G <- data.frame(R2G)
row.names(R2G) <- paste0(R2G$seqnames,":",R2G$start-1,"-",R2G$end) ##-1 in Peaks
head(R2G)

table(row.names(beds) %in% row.names(R2G)) ##15 genes not present

##adding gene names
##loading ENS ids and getting gene name associated with them
ens <- read.delim("/path/to/Files/gencdH38p13r37CR_genes.txt",header = FALSE)
colnames(ens) <- c("Gene.ID","Gene_name")
ens <- ens %>% separate(Gene.ID,c("Gene.ID"))
head(ens)
ens <- ens[!duplicated(ens$Gene.ID),]
row.names(ens) <- ens$Gene.ID
R2G$Gene.ID <- R2G$Genes <- NA

ids <- c("EOMES","CRX","Both")
GSEA <- list()
for(x in ids){
  peaks <- beds[beds$ID==x,"V4"]
  gene_sel <- c()
  for(y in peaks){
    gx <- unlist(R2G[y,"annotated_genes"])
    gene_names <- ens[gx,"Gene_name"]
    gene_sel <- c(gene_sel,gene_names)
  }
  GSEA[[x]] <- unique(gene_sel)
  rm(gene_sel,peaks)
}
rm(x)
GSEAx <- GSEA
GSEA[["EOMES"]] <- unique(c(GSEA$EOMES,GSEA$Both))
GSEA[["CRX"]] <- unique(c(GSEA$CRX,GSEA$Both))
GSEA[["Both"]] <- intersect(GSEA$EOMES,GSEA$CRX)
lengths(GSEA)
lengths(GSEAx)

save(GSEA, GSEAx,file = "TargetsofCRE_EOMES_CRX.RData")
