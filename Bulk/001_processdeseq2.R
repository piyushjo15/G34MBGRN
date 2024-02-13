##DESeq2 normalization of three data sets

suppressPackageStartupMessages({
    library(DESeq2)
})
##ICGC----
load("rawcounts_com.RData")
keep <- meta$Batch=="ICGC"
table(keep)
mdt.ICGC <- mdt[,keep]
ICGC <- meta[keep,]
rm(keep)
##DESeq2
dds.ICGC <- DESeqDataSetFromMatrix(countData = mdt.ICGC,
                              colData = ICGC,
                              design = ~ Subgroup)
dds.ICGC <- DESeq(dds.ICGC)
vsd.ICGC <- vst(dds.ICGC, blind=FALSE)
mdt.vsd.ICGC <- assay(vsd.ICGC)

tvg.ICGC <-row.names(mdt.vsd.ICGC)[order(rowVars(assay(vsd.ICGC)),decreasing=TRUE)]
##remove duplciated genes here
genes.id <- read.delim("geneidhg38r37.txt",row.names = 1)
genes.id <- genes.id[tvg.ICGC,]
keep <- duplicated(genes.id$Gene_name)
table(keep)
genes.id <- genes.id[!keep,]
tvg.ICGC <- genes.id$Gene_name
mdt.vsd.ICGC <- mdt.vsd.ICGC[row.names(genes.id),]
row.names(mdt.vsd.ICGC) <- genes.id$Gene_name

mdt.ICGC <- mdt.ICGC[row.names(genes.id),]
row.names(mdt.ICGC) <- genes.id$Gene_name

save(dds.ICGC, tvg.ICGC, mdt.ICGC, mdt.vsd.ICGC, ICGC, meta,file = "ICGC_DS2.RData")
rm(dds.ICGC, tvg.ICGC, mdt.ICGC, mdt.vsd.ICGC, keep)
##CR ----
keep <- meta$Batch%in% c("A1905","Eurofins")
table(keep)
mdt.CR <- mdt[,keep]
CR <- meta[keep,]
rm(keep)

mdt.CR <- mdt.CR[row.names(genes.id),]
row.names(mdt.CR) <- genes.id$Gene_name

##DESeq2
dds.CR <- DESeqDataSetFromMatrix(countData = mdt.CR,
                                   colData = CR,
                                   design = ~ Subgroup+Batch)
dds.CR <- DESeq(dds.CR)
vsd.CR <- vst(dds.CR, blind=FALSE)
mdt.vsd.CR <- assay(vsd.CR)

tvg.CR <-row.names(mdt.vsd.CR)[order(rowVars(assay(vsd.CR)),decreasing=TRUE)]
save(dds.CR, tvg.CR, mdt.CR, mdt.vsd.CR,CR, file = "CR_DS2.RData")
rm(dds.CR, tvg.CR, mdt.CR, mdt.vsd.CR,keep)
##MDT ----
keep <- meta$Batch%in% c("MDT")
table(keep)
mdt.MDT <- mdt[,keep]
MDT <- meta[keep,]

mdt.MDT <- mdt.MDT[row.names(genes.id),]
row.names(mdt.MDT) <- genes.id$Gene_name

##DESeq2
dds.MDT <- DESeqDataSetFromMatrix(countData = mdt.MDT,
                                 colData = MDT,
                                 design = ~ Subgroup)
dds.MDT <- DESeq(dds.MDT)
vsd.MDT <- vst(dds.MDT, blind=FALSE)
mdt.vsd.MDT <- assay(vsd.MDT)

tvg.MDT <-row.names(mdt.vsd.MDT)[order(rowVars(assay(vsd.MDT)),decreasing=TRUE)]
save(dds.MDT, tvg.MDT, mdt.MDT, mdt.vsd.MDT, MDT,file = "MDT_DS2.RData")
rm(dds.MDT, tvg.MDT, mdt.MDT, mdt.vsd.MDT,keep)
q()