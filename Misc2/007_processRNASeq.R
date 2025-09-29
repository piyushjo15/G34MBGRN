## This script process RNA-seq counts data from HDMB03/MB3W1 cell-line/PDX samples
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(BiocParallel)
  library(tidyverse)
})

setwd("MBRS/")
#for feqature counts----
## data avialble on https://doi.org/10.1101/2024.02.09.579680

mat <- read.delim("HDMB03_PDX_mat.txt")
genes <- row.names(mat)
grp ## metadata
table(colnames(mat) == row.names(grp))

dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = grp,
                              design = ~BATCH+Cond)
dds
# adding gene name info to DESeq metadata
mcols(dds) <- DataFrame(mcols(dds), genes)

keep <- rowSums(counts(dds)) >= 3
dds <- dds[keep,]

dds <- DESeq(dds, parallel = TRUE)

vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)
save(mat,vsd, grp,dds, file = "normcount_HDMB03_PDX.RData")

##log fold analysis-------------
load("normcount_HDMB03_PDX.RData")
resultsNames(dds)
res <- lfcShrink(dds, contrast=c("Cond","CRX_KO","CTRL"),type="ashr") 

# res <- lfcShrink(dds, contrast=c("Cond","EOMES_OE","CTRL"),type="ashr")
# 
# res <- lfcShrink(dds, contrast=c("Cond","MYC_KD","CTRL"),type="ashr")
# 
# res <- lfcShrink(dds, contrast=c("Cond","PAX6_OE","CTRL"),type="ashr") 

# res <- lfcShrink(dds, contrast=c("Cond","EOMES_OE_MYC_KD","MYC_KD"),type="ashr")

res2 <- data.frame(res)
res2$Gene <- row.names(res2)
res2 <- res2[order(res2$pvalue),]
head(res2)

##remove genes with NAN
res2 <- res2[!is.na(res2$padj),]
res2$LogFDR <-(-1) *log10(res2$padj)

res2[res2$LogFDR=="Inf","LogFDR"] <- max(res2$LogFDR[res2$LogFDR!="Inf"])
res2$significant <- with(res2, ifelse(LogFDR >= 1.3 & log2FoldChange > 0.5, "Up",
                                      ifelse(LogFDR >= 1.3 & log2FoldChange < -0.5, "Down", "No")))
##capping folder change
res2[res2$log2FoldChange>4,"log2FoldChange"] <-4
res2[res2$log2FoldChange<(-4),"log2FoldChange"] <-(-4)
res2[res2$LogFDR>10,"LogFDR"] <-10

table(res2$significant)
res2$significant2 <- res2$significant


##load gene-sets
load("Merged_GSEA.RDATA") ## axial gene-set

### comparing to Axial markers
kk <- res2$Gene[res2$significant %in% c("Up","Down")]
sel_pr <- GSEA$PR
sel_ubc <- unique(c(GSEA$Prec,GSEA$UBC))

## on signficant geens
sel_ubc <- sel_ubc[sel_ubc %in% kk]
sel_pr <- sel_pr[sel_pr %in% kk]

res2[res2$Gene %in% sel_pr,"significant2"] <- "G3"
res2[res2$Gene %in% sel_ubc,"significant2"] <- "G4"

tit <- "HDMB03 CRX_KO PDX PR vs UBC"

res2$significant2 <- factor(res2$significant2,levels = c("No","Down","Up","G3","G4"))
p <- ggplot(res2 %>% arrange(significant2), aes(x = log2FoldChange, y = LogFDR)) +
  geom_point(aes(color = significant2),size=2) +
  scale_color_manual(values = c("grey87","tan1","lightblue", "red4","darkgreen")) +
  geom_vline(xintercept = c(-.5, .5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(axis.title = element_text(size=10),
        axis.text= element_text(size=10),
        panel.grid.minor = element_blank(),
        legend.position = "None")+
  scale_x_continuous(limits = c(-4.1,4.1))+
    labs(title = tit,
       x = "Log2 Fold Change",
       y = "-Log10(FDR)",
       color = "Significance")
tit2 <- gsub(" ","_",tit)

pdf(paste0("Figures/",tit2,".pdf"),  width=3, height=3, pointsize = 10)
print(p)
dev.off()


##make a bar plot

df <-data.frame(Genes=c(sel_pr,sel_ubc),
                Group=c(rep("Group 3",length(sel_pr)),rep("Group 4",length(sel_ubc))),
                Dir="Up")
dwn <- res2$Gene[res2$significant %in% c("Down")]
df[df$Genes %in% dwn,"Dir" ] <- "Down"
df$Group <- factor(df$Group,levels = c("Group 3","Group 4"))
df$Dir <- factor(df$Dir,levels = c("Down","Up"))
vx <- c("tan1","lightblue")
names(vx) <- c("Down","Up")
p <- ggplot(df, aes(x =Group, fill=Dir)) +
  geom_bar(position = "fill")+
  scale_fill_manual(values =vx ) +
  theme_minimal() +
  theme(axis.text= element_text(size=10),
        axis.title = element_blank(),
        legend.position = "None")

pdf(paste0("Figures/",tit2,"_bar.pdf"),  width=2, height=3, pointsize = 10)
print(p)
dev.off()



#plotting PCA ----
#vsd <- vst(vsd, blind=FALSE)
# to change top HVG gene considered for PCA.
# check 500, 1000, 1500 (or more)
p <- plotPCA(vsd, intgroup=c("Cond"),ntop = 1500 )
d <- length(unique(grp$Cond))
p2 <- p + geom_point(aes(fill = grp$Cond), size=2)+
  geom_text_repel(aes(label=name), size=2)+ # comment this out if you don't want labels
  scale_color_manual(values =colorRampPalette(pal_simpsons()(14))(d))+ # to change point color
  theme_classic()+
  theme(axis.line = element_line(colour = 'black',linewidth=.5),
        axis.ticks = element_line(colour = 'black',linewidth=0.5),
        axis.text.x = element_text(face="bold",colour = "black",size=8),
        axis.text.y = element_text(face="bold",colour = "black",size=8),
        axis.title.y=element_text(face="bold",colour="black",size = 10),
        axis.title.x=element_text(face="bold",colour="black",size = 10),
        legend.position = "None")
#to save as plot tiff file
png("Figures/HDMB03_PDX.png", units="in", width=8, height=8, res=300)
print(p2)
dev.off()

plotPCA(vsd, "Batch")

assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$REP)

plotPCA(vsd, "Batch")

``# plotting gene exp boxplot-------
load("genematrix_HDMB03_PDX.RData")

mat <- mat[,row.names(grp)]
table(colnames(mat) == row.names(grp))
##for counts
keep <- duplicated(genes$name)
genes <- genes[!keep,]
mat <- as.matrix(mat[!keep,])
row.names(mat) <- genes$name
mat[1:4,1:4]


plotbx <- function(x){
  #pd <- grp %>% mutate(Gene=mat[x,]) #VSD
  #pd <- grp %>% mutate(Gene=log2(mat[x,]+1)) #counts
  pd <- grp %>% mutate(Gene=(mat[x,]+1))
  pd$Cond <- factor(pd$Cond,levels = sort(unique(pd$Cond)))
  p <- ggplot(pd,aes(x=Cond,y=Gene ))+
    geom_boxplot()+
    geom_text(label=row.names(pd),size=4)+
    geom_point()+
    theme_light()+
    theme(axis.line = element_line(colour = 'black',linewidth = 0.5),
          axis.ticks = element_line(colour = 'black',linewidth=0.5),
          axis.text.x = element_text(face="bold",colour = "black",size=8, angle = 90,vjust=0.5, hjust=1),
          axis.text.y = element_text(face="bold",colour = "black",size=8),
          axis.title = element_text(face="bold",colour = "black",size=10),
          legend.position = "None")+
    labs(y="Exp. Counts:", x="Experimental conditions")+
    ggtitle(x)
  return(p)
}

x <- "LMX1A"
png(paste0("Figures/HDMB03_PDX_",x,".png"),units="in", width = 6, height = 5, res=300)
plotbx(x)

dev.off()

x <- "PVT1"

plotbx(x)

df <- mat[c("CRX","EOMES","LMX1A","PAX6","LHX1","TRPC3","NRL","NR2E3","SAG"),row.names(grp)[grp$Cond=="PAX6_OE"]]
df <- mat[c("CRX","EOMES","PAX6"),row.names(grp)[grp$Cond=="PAX6_OE"]]
df <- log2(df+1)
df2 <- data.frame(t(df))
ggplot(df2,aes(x=SAG,y=TRPC3))+geom_point()+geom_text_repel(label=row.names(df2))+theme_bw()
pairs(df2,pch=19,lower.panel = NULL)
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19)
}
# Create the plots
pdf(paste0("Figures/HDMB03_PDX_COR.pdf"), width = 3, height =3, pointsize = 10)
pairs(as.matrix(df2), 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
dev.off()

load("genematrix_MB3W1_PDX.RData")


grp <- grp[grp$Type=="MB3W1",]
plotbx2 <- function(sel){
  genes <- cand
  del <- apply(mat[genes,]+1, 1, log10)
  pd <- cbind(grp[,c("Batch","Cond")],del)
  pd <- pd[pd$Cond %in% c("CTRL",sel),]
  pd$Batch <- NULL
  pd2 <- reshape2::melt(pd)
  rm(pd,del)
  pd2$variable <- factor(pd2$variable,levels = genes)
  p <- ggplot(pd2,aes(x=variable,y=value, fill=Cond, color=Cond))+
    geom_boxplot(position=position_dodge(width = 1))+
    geom_point(position=position_dodge(width = 1))+
    scale_color_manual(values = c("grey80","black"))+
    scale_fill_manual(values = c("grey80","black"))+
    theme_classic()+
    theme(axis.line = element_line(colour = 'black',linewidth = 0.5),
          axis.ticks = element_line(colour = 'black',linewidth=0.5),
          axis.text.x = element_text(face="bold",colour = "black",size=8, angle = 90,vjust=0.5, hjust=1),
          axis.text.y = element_text(face="bold",colour = "black",size=8),
          axis.title = element_text(face="bold",colour = "black",size=10))+
    labs(y="Exp. Counts:")
  return(p)
}

cand <- c("MYC","FOXN4","E2F7","PBX3",
          "CRX","SAG","NR2E3","NRL","RXRG",
          "EOMES","LMX1A","SNCAIP","RPH3A","GRM8",
          "PAX6")

sel <- "PAX6_OE"
pdf(paste0("MBRS/HDMB03CL_",sel,".pdf"), width = 6, height = 3.5, pointsize = 10)
plotbx2(sel)
dev.off()
#####
setwd("~/MBsnANA/UBC_TF_ana/MBRS/")

load("~/MBnewana/PSBHDMB03.RData")
mat2 <- logcounts(del)
load("normcount_HDMB03_PDX.RData")
grp <- grp[grp$Cond %in% c("CTRL","PAX6_OE"),]
mat <- mat[,row.names(grp)]
# Get the names of the top 20 most variable genes
hvg <- order(rowVars(mat), decreasing = TRUE)[1:4000]
hvg <- row.names(mat)[hvg]

hvg <- intersect(hvg[1:2000],row.names(mat2))
mat1x <- t(apply(mat, 1, scale))

colnames(mat1x) <- colnames(mat)

mat2x <- t(apply(mat2, 1, scale))
colnames(mat2x) <- colnames(mat2)


corx <- cor(mat1x[hvg,],mat2x[hvg,])

corrplot::corrplot(corx,method = "number",tl.cex = .4)


