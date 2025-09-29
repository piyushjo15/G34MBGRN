## GSEA analysis for Gene-sets

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(DESeq2)
  library(fgsea)
})

setwd("~MBRS/")
##load gene-sets
load("Merged_GSEA.RDATA")

GSEAx <- list()
GSEAx[["PR"]] <- GSEA$PR
GSEAx[["Prec"]] <- GSEA$Prec
GSEAx[["UBC"]] <- GSEA$UBC
GSEA <- GSEAx

lengths(GSEA)
rm(GSEAx)

## 1. load normalized counts ----------
## HDMB03 PDX -------------
load("normcount_HDMB03_PDX.RData")
resultsNames(dds)

res <- lfcShrink(dds, contrast=c("Cond","CRX_KO","CTRL"),type = "ashr") 
tit <- "HDMB03 CRX_KO PDX"

res <- lfcShrink(dds, contrast=c("Cond","EOMES_OE","CTRL"),type = "ashr") 
tit <- "HDMB03 EOMES_OE PDX"

res <- lfcShrink(dds, contrast=c("Cond","PAX6_OE","CTRL"),type = "ashr") 
tit <- "HDMB03 PAX6_OE PDX"

## HDMB03 CL -------------
load("normcount_HDMB03_CL.RData")
resultsNames(dds)

res <- lfcShrink(dds, contrast=c("Cond","CTRL","WT"),type = "ashr") 
tit <- "HDMB03 CTRL CL"

res <- lfcShrink(dds, contrast=c("Cond","MYC_KD","CTRL"),type = "ashr") 
tit <- "HDMB03 MYC_KD CL"

res <- lfcShrink(dds, contrast=c("Cond","EOMES_OE","CTRL"),type = "ashr") 
tit <- "HDMB03 EOMES_OE CL"

res <- lfcShrink(dds, contrast=c("Cond","PAX6_OE","CTRL"),type = "ashr") 
tit <- "HDMB03 PAX6_OE CL"

res <- lfcShrink(dds, contrast=c("Cond","EOMES_OE_MYC_KD","MYC_KD"),type = "ashr") 
tit <- "HDMB03 EOMES_OE MYC_KD CL"

## MB3W1 CL -------------
load("normcount_MB3W1_Cell.RData")
resultsNames(dds)
res <- lfcShrink(dds, contrast=c("Cond","MYC_KD","CTRL"),type="ashr") 
tit <- "MB3W1 MYC_KD CL"

res <- lfcShrink(dds, contrast=c("Cond","EOMES_OE","CTRL"),type="ashr")
tit <- "MB3W1 EOMES_OE CL"

res <- lfcShrink(dds, contrast=c("Cond","EOMES_OE_MYC_KD","MYC_KD"),type="ashr") 
tit <- "MB3W1 EOMES_OE MYC KD CL"

res <- lfcShrink(dds, contrast=c("Cond","PAX6_OE","CTRL"),type="ashr") 
tit <- "MB3W1 PAX6_OE CL"


## MB3W1 PDX -------------
load("normcount_MB3W1_PDXv2.RData")
resultsNames(dds)
res <- lfcShrink(dds, contrast=c("Cond","MYC_KD","CTRL"),type="ashr") 
tit <- "MB3W1 MYC_KD PDX"

res <- lfcShrink(dds, contrast=c("Cond","EOMES_OE","CTRL"),type="ashr")
tit <- "MB3W1 EOMES_OE PDX"

res <- lfcShrink(dds, contrast=c("Cond","PAX6_OE","CTRL"),type="ashr") 
tit <- "MB3W1 PAX6_OE PDX"


## 2. processing for GSEA ----------------
### rank the gene list by -log10 fdr
res2 <- data.frame(res)
res2$Gene <- row.names(res2)
#head(res2)
##remove genes with NAN
res2 <- res2[!is.na(res2$pvalue),]
res2$LogFDR <-(-1) *log10(res2$pvalue)
res2$Sign <- 1
res2[res2$log2FoldChange<0,"Sign"] <- -1
res2$LogFDR2 <- res2$Sign * res2$LogFDR
res2 <- res2[order(res2$LogFDR2,decreasing = TRUE),]

genes <- res2$LogFDR2
names(genes) <- res2$Gene
rm(res,res2)

set.seed(1)
ff <- fgseaMultilevel(
  pathways = GSEA,
  stats = genes,
  scoreType='std',
  nproc=6
)

ff

values <- read.delim("Files/axiscol.txt",row.names = 1)
cls <- sort(names(GSEA))
valuesx <- values[cls,1]


pd <- c()
for(x in cls){
  pd_del <- plotEnrichmentData(GSEA[[x]], genes)
  pd_del  <- pd_del$curve
  pd_del$Class <- x
  pd <- rbind(pd,pd_del)
  rm(pd_del)
}


p <- ggplot(pd,aes(x=rank,y=ES,color=Class))+
  geom_line(linewidth=1)+
  scale_color_manual(values = valuesx)+
  theme_minimal()+
  geom_hline(yintercept = 0,color="grey34",linetype="dashed",linewidth=2)+
  theme(axis.text = element_text(size=10),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 2),
        legend.position = "None")+
  scale_x_continuous(label=function(x) {
    parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
  }
  )

tit2 <- gsub(" ","_",tit)


pdf(paste0("Figures/",tit2,"_GSEA.pdf"),  width=4, height=3, pointsize = 10)
print(p)
dev.off()

