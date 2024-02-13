## This script uses AUC score obtained from Bulk RNA-Seq samples and creates
## TF-GRN enrichment heatmap per subtype and obtain TSNE map based on TF-GRN scores

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

library(scater)
library(RColorBrewer)
load("pdMB_G34aucGRN_AUC.RData")
ms2 <- t(assay(ms))
sel <- read.delim("GRNs.txt")
sel <- sel$GRNs

# ##WGCNA
# load("pdMB_WGCNA_AUC.RData")
# ms2 <- scr
# sel <- colnames(scr)
#remove unannotated samples------
keep1 <- is.na(meta$Subtype_M)
keep3 <- meta$Subgroup %in% c("SHH","WNT")
keep2 <- meta$Subtype_M %in% c("MYO_NA","MB_MYO","RB")
keep <- keep1 | keep2 | keep3
table(keep)
ms1 <- ms2[!keep,]
metax <- meta[!keep,]
rm(keep1,keep2,keep)

table(metax$Subtype_M)
table(metax$Subgroup_M)
pd <- metax[,c(4,3,6)]

ms1x <- ms1[,sel]
##calculate zscores separately-----
keep <- pd$Batch=="ICGC"
table(keep)
ms1icgc <- ms1x[keep,]
ms1icgc2 <- apply(ms1icgc, 2, zscore_c)

keep <- pd$Batch=="MDT"
table(keep)
ms1mdt <- ms1x[keep,]
ms1mdt2 <- apply(ms1mdt, 2, zscore_c)

keep <- pd$Batch %in% c("A1905","Eurofins")
table(keep)
ms1cr <- ms1x[keep,]
ms1cr2 <- apply(ms1cr, 2, zscore_c)

ms2 <- rbind(ms1icgc2,ms1cr2,ms1mdt2)
save(ms1x, ms2, pd, file = "modscorecomzAUC.RData")

# #WCGNA
# save(ms1x, ms2, pd, file = "WCGNA_com_zAUC.RData")


#ann color-----
load("modscorecomzAUC.RData")

###-------
corx <- cor(ms1x, ms1x)
myColor <- colorRampPalette(c("blue","white","red4"))(100)

myBreaks <- c(seq(min(corx), 0, length.out=ceiling(100/2) + 1), 
              seq(max(corx)/100, max(corx), length.out=floor(100/2)))
hp <- pheatmap::pheatmap(corx,clustering_method = "ward.D2",
                         clustering_distance_rows = "correlation",clustering_distance_cols =  "correlation",
                         color=myColor, breaks=myBreaks,
                         border_color = NA, 
                         fontsize = 8)
tiff(filename = "G34MB_cor_GEP_D.tiff", width = 18, height = 18, units = "in", res = 300)
print(hp)
dev.off()

##enr-----

ds <- read.delim("dataset.txt",row.names = 1)
dsc <- ds$Color
names(dsc) <- row.names(ds)

gs <- read.csv("MBcolorG34.csv",row.names = 1)
gsc <- gs$Color
names(gsc) <- row.names(gs)

subs <- read.delim("Subtype_G34.txt", row.names = 1)
sc <- subs$Color
names(sc) <- row.names(subs)

ann_col <- list(Batch=dsc,
                Subgroup_M=gsc,
                Subtype_M=sc)

pd$Batch <- factor(pd$Batch, levels = row.names(ds))
pd$Subgroup_M <- factor(pd$Subgroup_M, levels = row.names(gs))
pd$Subtype_M <- factor(pd$Subtype_M, levels = row.names(subs))
##-----
## plotting TF-GRN heatmap in tumor samples per subtype
subsx <- row.names(subs)
for(x in subsx){
  keeps <- pd$Subtype_M==x
  table(keeps)
  pds <- pd[keeps,]
  ms2s <- ms2[keeps,]
  h <- round(dim(pds)[1]*.05)
  hp <- pheatmap::pheatmap(ms2s[,sel], clustering_method = "ward.D2",
                           color=myColor, breaks=myBreaks,clustering_distance_rows = "correlation",
                           cluster_cols = FALSE,
                           fontsize = 7,legend = FALSE,annotation_legend = FALSE,
                           border_color = NA, show_rownames = FALSE,show_colnames = FALSE,
                           treeheight_row =  FALSE,annotation_names_row = FALSE,
                           annotation_row = pd, annotation_colors = ann_col)
  
  tiff(filename = paste0(x,"_enr_E.tiff"), width = 8, height = h, units = "in", res = 300)
  print(hp)
  dev.off()
  # dela <- hp$tree_row
  # delb <- dela$labels
  # order_g1 <- delb[dela$order]
  # write(order_g1, paste0(x,"_og.txt"))
  #rm(dela,delb,order_g1,hp, keeps, pds, ms2s,h)
  rm(hp, keeps, pds, ms2s,h)
}
rm(x)


##selected samples
sel <- paste0("ICGC_MB",c(129,292,26))
ms2x <- ms2[sel,]
hp <- pheatmap::pheatmap(ms2x, 
                         color=myColor, breaks=myBreaks,
                         cluster_cols = FALSE,cluster_rows = FALSE,
                         fontsize = 7,
                         border_color = NA, 
                         #show_rownames = FALSE,
                         show_colnames = FALSE,
                         treeheight_row =  FALSE)

tiff(filename = "G34MB_sel_VII.tiff", width = 4, height = .5, units = "in", res = 300)
print(hp)
dev.off()

##TSNE ----
plot.data <- pd

set.seed(345)
tsne <- Rtsne::Rtsne(ms2, pca=FALSE,perplexity=30, max_iter=5000,pca_center=FALSE,
                     normalize=FALSE, theta=0.0)
plot.data$TSNE1 <- tsne$Y[,1]
plot.data$TSNE2 <- tsne$Y[,2]

tiff( "TSNE.tiff", width = 4, height = 4, units = "in", res = 300)
ggplot(plot.data, aes(TSNE1, TSNE2, color=Subtype_M))+
  geom_point(size=2)+
  scale_color_manual(values=sc)+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title=element_blank(),
        legend.position = "None")
dev.off()


ggplot(plot.data, aes(TSNE1, TSNE2, color=Batch))+
  geom_point(size=2)+
  scale_color_manual(values=dsc)+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title=element_blank(),
        legend.position = "None")

save(plot.data, file = "pdTSNG34MB.RData")

## plotting individual gene-set scores on TSNE ---------------
load("modscorecomzAUC.RData")
#load("WCGNA_com_zAUC.RData")
load("pdTSNG34MB.RData")


plotMS <- function(x){
  MS=scale(range_x(ms1x[,x]))
  MS[MS>2] <- 2
  MS[MS<(-2)] <- (-2)
  pdel <- cbind(plot.data,MS=MS )
  pr <- ggplot(pdel %>% arrange(MS), aes(x=TSNE1,y=TSNE2, color=MS))+
    geom_point()+
    scale_color_gradientn(colours = colorRampPalette(c("deepskyblue4","beige","magenta3"))(100),
                          breaks=c(-2,0,2),
                          limits=c(-2, 2))+
    theme_classic()+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title=element_blank(),
          legend.position = "None")
  return(pr)
  
}
x <- "ALX1"
tiff(paste0(x,"_TSNE.tiff") , width = 4, height = 4, units = "in", res = 300)

plotMS(x)
dev.off()

#WGCNA
sel2 <- c("green","black",
          "midnightblue","brown")
#### plotting scaled gene expression score on TSNE-------
##zscore approach----
load("ICGC_DS2.RData")
rm(ICGC,dds.ICGC,mdt.ICGC,meta,tvg.ICGC)
load("MDT_DS2.RData")
rm(MDT,dds.MDT,mdt.MDT, tvg.MDT)
load("CR_DS2.RData")
rm(CR,dds.CR,mdt.CR, tvg.CR)
mdt.vsd.ICGC2 <- apply(t(mdt.vsd.ICGC), 2, zscore_c)
mdt.vsd.CR2 <- apply(t(mdt.vsd.CR), 2, zscore_c)
mdt.vsd.MDT2 <- apply(t(mdt.vsd.MDT), 2, zscore_c)
rm(mdt.vsd.ICGC,mdt.vsd.CR,mdt.vsd.MDT)

vsd.all <- rbind(mdt.vsd.ICGC2,mdt.vsd.CR2,mdt.vsd.MDT2)
save(vsd.all, file="zscoreVSDall.RData")
## 
load("zscoreVSDall.RData")
load("pdTSNG34MB.RData")

vsd.all <- vsd.all[row.names(plot.data),]

plotgenes <- function(x){
  gx <- vsd.all[,x]
  gx[gx>2] <- 2
  gx[gx<(-2)] <- (-2)
  pdel <- cbind(plot.data,Gene=gx)
  pr <- ggplot(pdel %>% arrange(Gene), aes(x=TSNE1, y=TSNE2, color=Gene))+
    geom_point()+
    scale_color_gradientn(
      colours = colorRampPalette(c("skyblue","grey93","orangered4"))(100),
      breaks=c(-2,0,2),
      limits=c(-2, 2))+
    theme_classic()+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title=element_blank(),
          legend.position = 'None')
  return(pr)
}

x <- "PRDM6"
tiff(paste0(x,"_TSNE_genesz.tiff") , width = 4, height = 4, units = "in", res = 300)
plotgenes(x)
dev.off()
