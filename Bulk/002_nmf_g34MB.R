## This script is to obtain meta gene score for bulk data to show continumm
library(reticulate)

Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")

#function for  return nearest neighbour less than radius and max neighbors
suppressPackageStartupMessages({
  library(Matrix)
  library(uwot)
  library(DESeq2)
  
})
##load data-----
load("ICGC_DS2.RData")
load("CR_DS2.RData")
load("MDT_DS2.RData")
rm(list=ls(pattern = "^dds*"))
rm(mdt.ICGC, mdt.CR, mdt.MDT, meta)
genes <- rownames(mdt.vsd.ICGC)
mdt <- cbind(mdt.vsd.ICGC,mdt.vsd.CR,mdt.vsd.MDT)
meta <- rbind(ICGC,CR,MDT)
rm(CR, ICGC, MDT, mdt.vsd.CR, mdt.vsd.MDT, mdt.vsd.ICGC)
##remove the non important samples----
table(row.names(meta)==colnames(mdt)) ##alll check
keep1 <- is.na(meta$Subtype_M)
keep3 <- meta$Subgroup %in% c("SHH","WNT")
keep2 <- meta$Subtype_M %in% c("MYO_NA","MB_MYO","RB")
keep <- keep1 | keep2 | keep3
table(keep)
metax <- meta[!keep,] ##703 G3/4 samples
rm(keep1,keep2,keep, meta)

###identify common tvgs 
tvg <- intersect(tvg.CR[1:5000],intersect(tvg.ICGC[1:5000],tvg.MDT[1:5000]))
rm(tvg.CR, tvg.ICGC, tvg.MDT)
#subset
mdt <- mdt[tvg,row.names(metax)]
table(row.names(metax)==colnames(mdt)) ##alll check
dim(mdt)
table(metax$Batch)
###preparing for NMF-----
##using MDT samples for reference in this approach
keep <- metax$Batch=="MDT"
table(keep)
mdt1 <- t(mdt[,keep])
mdt2 <- t(mdt[,!keep])

##python run ----
repl_python()
import sklearn.decomposition as sk
import numpy as np
import random

model1 = sk.NMF(n_components=2,init="nndsvd", 
               max_iter=100000, random_state=0, tol=1e-5) 
random.seed(456)
model1.fit(r.mdt1)
wna1 = model1.transform(r.mdt1)
wnb1 = model1.transform(r.mdt2)
h1=model1.components_

model2 = sk.NMF(n_components=8,init="nndsvd", 
                max_iter=100000, random_state=0, tol=1e-5) 
random.seed(456)
model2.fit(r.mdt1)
wna2 = model2.transform(r.mdt1)
wnb2 = model2.transform(r.mdt2)
h2=model2.components_

exit

rda_2 <- py$wna1
rdb_2 <- py$wnb1
h1 <- py$h1

rda_8 <- py$wna2
rdb_8 <- py$wnb2
h2 <- py$h2

row.names(h1) <-colnames(rda_2) <- colnames(rdb_2) <-paste0("NMF_",seq(dim(rda_2)[2]))
row.names(h2) <-colnames(rda_8) <- colnames(rdb_8) <-paste0("NMF_",seq(dim(rda_8)[2]))
row.names(rda_2) <-row.names(rda_8) <- row.names(mdt1)
row.names(rdb_2) <- row.names(rdb_8) <- row.names(mdt2)

colnames(h1) <-colnames(h2) <- tvg
rd_2 <- rbind(rda_2,rdb_2)
rd_2 <- rd_2[row.names(metax),]

rd_8 <- rbind(rda_8,rdb_8)
rd_8 <- rd_8[row.names(metax),]
##umap
set.seed(134)
del1 <- umap(rd_8, n_components = 2, min_dist = 0.8)

metax$nfUMAP1 <- del1[,1]
metax$nfUMAP2 <- del1[,2]
save(rd_8,rd_2, h1,h2, metax, mdt, tvg, file ="NMF_G34bulkRS.RData" )
q()
##post NMF ----
load("NMF_G34bulkRS.RData")
##plots------------
##batch

d <- sort(unique(metax$Batch))
ds <- read.delim("Files/dataset.txt",row.names = 1)
dsc <- ds[d,1]
p <- ggplot(metax, aes(x=nfUMAP1, nfUMAP2, color=Batch))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values = dsc)+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
tiff( "UMAP_Batch_NMF.tiff", width = 6, height = 6, units = "in", res = 300)
print(p)
dev.off()
rm(p)
##subgrpup
groupx <- c("goldenrod1","darkgreen")
p <- ggplot(metax, aes(x=nfUMAP1, nfUMAP2, color=Subgroup_M))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values = groupx)+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
tiff( "UMAP_Group_NMF.tiff", width = 6, height = 6, units = "in", res = 300)
print(p)
dev.off()
rm(p)
##subtype
d <- sort(unique(metax$Subtype_M))
subtype <- read.delim("Files/Subtype_G34.txt",row.names = 1)
subtypex <- subtype[d,1]

p <- ggplot(metax, aes(x=nfUMAP1, nfUMAP2, color=Subtype_M))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values = subtypex)+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

tiff( "UMAP_Subtype_NMF.tiff", width = 6, height = 6, units = "in", res = 300)
print(p)
dev.off()
rm(p)
###finding markers -----
marks <- scran::findMarkers(t(rd_8), groups=metax$Subgroup_M, direction="up",
                     test.type="wilcox",pval.type="some")
# subs <- names(marks)
# sig <- c()
# for(x in subs){
#   del <- marks[[x]]
#   del <- row.names(del)[1:2]
#   sig <- rbind(sig,del)
#   rm(del)
# }
# row.names(sig) <- subs
# sigx <- unique(c(sig))

sigsx <- row.names(marks$G3)

####################
range_01 <- function(x){
  mx <- quantile(x,0.99)
  mn <- quantile(x,.1)
  x[x>mx] <- mx
  x[x<mn] <- mn
  a <- 1/(max(x)-min(x))
  b <- 1-(a*max(x))
  d <- round((a*x)+b, digits = 3)
  return(d)
  
}
rdx8 <- apply(rd_8, 2, range_01)
rdx2 <- apply(rd_2, 2, range_01)
sigsx2 <- paste0("NMF_",c(3,8,2,7,5,6,4,1))
###arranging samples by scale
methscr <- read.delim("G34_methids.annotated.txt",row.names = 1)
head(methscr)
#load plot data
load("mnpscoresG34.RData")
plot.data$ID <- row.names(plot.data)
row.names(plot.data) <- plot.data$MethID
table(row.names(plot.data) %in% row.names(methscr))
methscr <- methscr[row.names(plot.data),]
row.names(methscr) <- plot.data$ID
table(row.names(metax) %in% row.names(methscr))
methscr <- methscr[row.names(metax),]
metax$G3_score <- methscr$MB..G3
metax$G4_score <- methscr$MB..G4

#metax$Scr <- rowMeans(rdx[,sigsx[1:3]]) -  rowMeans(rdx[,sigsx[6:8]])
metax$Scr <- rdx2[,2] -  rdx2[,1] ##G3(2)- G4(1) Metagene 
metax$Scrx <- methscr$MB..G3 -  methscr$MB..G4 ##G3- G4 Metagene 

##plotting metagene and methylations core-----
p1 <- ggplot(metax, aes(x=Scr, y=G3_score, color=Subgroup ))+
  geom_point()+
  scale_color_manual(values = groupx)+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black',size=.8),
        axis.ticks = element_line(colour = 'black',size=.8),
        axis.text.x = element_text(face="bold",colour = "black",size=10),
        axis.text.y = element_text(face="bold",colour = "black",size=10),
        axis.title=element_blank(),
        legend.position = "None")+
  scale_y_continuous(breaks = c(0,0.5,1))
  

tiff( "G3_meth_G34Metagene.tiff", width = 3, height = 2, units = "in", res = 300)
print(p1)
dev.off()

p1 <- ggplot(metax, aes(x=Scr, y=G4_score, color=Subgroup ))+
  geom_point()+
  scale_color_manual(values = groupx)+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black',size=.8),
        axis.ticks = element_line(colour = 'black',size=.8),
        axis.text.x = element_text(face="bold",colour = "black",size=10),
        axis.text.y = element_text(face="bold",colour = "black",size=10),
        axis.title=element_blank(),
        legend.position = "None")+
  scale_y_continuous(breaks = c(0,0.5,1))


tiff("G4_meth_G34Metagene.tiff", width = 3, height = 2, units = "in", res = 300)
print(p1)
dev.off()
####------
ord <- order(metax$Scr, decreasing = TRUE)

metaxd <- metax[ord,]
# metaxd <- metax[order(rdx[,1],rdx[,4],rdx[,6],rdx[,8]),]

ann <- metaxd[,c("Subtype_M","Subgroup_M","Batch")]

names(subtypex) <- row.names(subtype)
names(groupx) <- c("G3","G4")
names(dsc) <- row.names(ds)

ann_col <- list(Batch=dsc,
                Subgroup_M=groupx,
                Subtype_M=subtypex)


mc <- colorRampPalette(c("white",brewer.pal(9,"OrRd")))(100)
mb <- ((seq(1:101)-1)/100)
hp <- pheatmap::pheatmap(t(rdx2[row.names(metaxd),]), 
                   annotation_col = ann,annotation_colors =ann_col,  
                   show_colnames =  FALSE, cluster_rows = FALSE, 
                   cluster_col = FALSE,
                   breaks = mb, color = mc,
                   annotation_legend = FALSE, legend = FALSE)
tiff( "Signautre_NMF_G34.tiff", width = 12, height = 6, units = "in", res = 300)
print(hp)
dev.off()
# hp <- pheatmap::pheatmap(t(methscr[row.names(metaxd),c("MB..G3","MB..G4")]), 
#                          annotation_col = ann,annotation_colors =ann_col,  
#                          show_colnames =  FALSE, cluster_rows = FALSE, 
#                          cluster_col = FALSE,
#                          breaks = mb, color = mc,
#                          annotation_legend = FALSE, legend = FALSE)
# tiff( "Signautre_NMF_G34.tiff", width = 12, height = 6, units = "in", res = 300)
# print(hp)
# dev.off()
hp <- pheatmap::pheatmap(t(rdx8[row.names(metaxd),sigsx2]), 
                         annotation_col = ann,annotation_colors =ann_col,  
                         show_colnames =  FALSE, cluster_rows = FALSE, 
                         cluster_col = FALSE,
                         breaks = mb, color = mc,
                         annotation_legend = FALSE, legend = FALSE)
tiff( "Signature_NMF_subtype.tiff", width = 12, height = 6, units = "in", res = 300)
print(hp)
dev.off()
##average per subtype
subs <- row.names(subtype)
msd2 <- c()
for(x in subs){
  keep <- metax$Subtype_M==x
  msd <- rd_8[keep,]
  del <- colMeans(msd)
  msd2 <- rbind(msd2,del)
  rm(keep, msd, del)
  
}
rm(x)
row.names(msd2) <- subs
msd3 <- apply(msd2, 2, range_01)
hp <- pheatmap::pheatmap(msd3[,sigsx2], 
                         cluster_col = FALSE,cluster_rows = FALSE,
                         breaks = mb, 
                         border_color = "black",
                         color = mc)
tiff( "Signautre_NMF2.tiff", width = 6, height = 4, units = "in", res = 300)
print(hp)
dev.off()

sigsx2 <- paste0("NMF_",c(3,8,2,7,5,6,4,1))
plotNMF <- function(x){
  pd <- cbind(metax[,c("Subgroup_M","Subtype_M")], NMF=rdx[,x])
  p <- ggplot(pd, aes(x=Subtype_M, y=NMF, fill=Subtype_M))+
    geom_boxplot()+
    scale_fill_manual(values =subtypex)+
    theme_classic()+
    theme(axis.line = element_line(colour = 'black',linewidth = .6),
        axis.ticks = element_line(colour = 'black',linewidth=.6),
        axis.text.x = element_text(face="bold",colour = "black",size=10, angle = 90, vjust = 0.5,hjust=1),
        axis.text.y = element_text(face="bold",colour = "black",size=10),
        axis.title=element_blank(),
        legend.position = "None")+
    scale_y_continuous(breaks = c(0,0.5,1))
  return(p)
}
for(i in 1:dim(rd_8)[2]){
  tiff(paste0("Signautre_NMF_",i,"_bp.tiff"), width = 4, height = 3, units = "in", res = 300)
  print(plotNMF(i))
  dev.off()
}

