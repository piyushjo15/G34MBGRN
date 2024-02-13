
suppressPackageStartupMessages({
  library(DESeq2)
})
##ICGC----
load("rawcounts_com.RData")
keep <- is.na(meta$Subtype_M)
meta <- meta[!keep,]
keep <- meta$Subgroup_M %in% c("G3","G4")
meta <- meta[keep,]
meta <- meta[!(meta$Subtype_M=="MYO_NA"),]
keep <- meta$Batch %in% c("Eurofins","A1905")
meta[keep,"Batch"] <- "Newcastle"
table(meta$Batch)
mdt <- mdt[,row.names(meta)]

genes.id <- read.delim("geneidhg38r37.txt",row.names = 1)
rn <- row.names(mdt)
rn <- rn[rn %in% row.names(genes.id)]
genes.id <- genes.id[rn,]
keep <- duplicated(genes.id$Gene_name)
table(keep)
genes.id <- genes.id[!keep,]
mdt <- mdt[!keep,]
row.names(mdt) <- genes.id$Gene_name

##DESeq2
dds <- DESeqDataSetFromMatrix(countData = mdt,
                                   colData = meta,
                                   design = ~ Subtype_M+Batch)
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,] 
dds <- DESeq(dds)
resultsNames(dds)
save(dds, file="intersubtypecomp.RData")
sams <- c("G34_I","G34_II","G34_III","G34_IV","G34_V","G34_VI","G34_VIII")
SVIIcomp <- list()
for(x in sams){
  SVIIcomp[[x]] <- lfcShrink(dds, contrast=c("Subtype_M","G34_VII",x), type="ashr", lfcThreshold=2)
}
sams <- c("G34_II","G34_III","G34_IV","G34_V","G34_VI","G34_VII","G34_VIII")
SIcomp <- list()
for(x in sams){
  SIcomp[[x]] <- lfcShrink(dds, contrast=c("Subtype_M","G34_I",x), type="ashr", lfcThreshold=2)
}
save(SVIIcomp,SIcomp, file = "SI_SVII_DEG.RData")

sams <- names(SIcomp)
SI <- data.frame(NULL)
for(x in sams){
  del <- data.frame(SIcomp[[x]])
  del <- del[!(is.na(del$padj)),]
  del <- del[(del$log2FoldChange>2 & del$padj<0.001),]
  del <- del[order(del$padj),]
  del <- del[1:250,]
  del$Gene <- row.names(del)
  del["PAX6",]
  del <- del[,c("Gene","log2FoldChange","padj")]
  row.names(del) <- seq(250)
  del$vsSubtype <- x
  SI <- rbind(SI,del)
  rm(del)
  
}
sams <- names(SVIIcomp)
SVII <- data.frame(NULL)
for(x in sams){
  del <- data.frame(SVIIcomp[[x]])
  del <- del[!(is.na(del$padj)),]
  del <- del[(del$log2FoldChange>2 & del$padj<0.001),]
  del <- del[order(del$padj),]
  del <- del[1:250,]
  del$Gene <- row.names(del)
  del <- del[,c("Gene","log2FoldChange","padj")]
  row.names(del) <- seq(250)
  del$vsSubtype <- x
  SVII <- rbind(SVII,del)
  rm(del)
}
rm(x)

write.table(SI, file = "SI_top250.txt",sep = "\t",quote = FALSE)
write.table(SVII, file = "SVII_top250.txt",sep = "\t",quote = FALSE)

genes <- c("PAX6","CRX","NRL","OTX2","LMX1A","EOMES")

sams <- names(SIcomp)
SI <- data.frame(NULL)
for(x in sams){
  del <- data.frame(SIcomp[[x]])
  del <- del[genes,]
  del$Gene <- row.names(del)
  del <- del[,c("Gene","log2FoldChange","padj")]
  row.names(del) <- seq(6)
  del$vsSubtype <- x
  SI <- rbind(SI,del)
  rm(del)
  
}
sams <- names(SVIIcomp)
SVII <- data.frame(NULL)
for(x in sams){
  del <- data.frame(SVIIcomp[[x]])
  del <- del[genes,]
  del$Gene <- row.names(del)
  del <- del[,c("Gene","log2FoldChange","padj")]
  row.names(del) <- seq(6)
  del$vsSubtype <- x
  SVII <- rbind(SVII,del)
  rm(del)
}
rm(x)
write.table(SI, file = "SI_sel.txt",sep = "\t",quote = FALSE)
write.table(SVII, file = "SVII_sel.txt",sep = "\t",quote = FALSE)
