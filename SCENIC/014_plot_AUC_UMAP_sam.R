### This script is for  is for getting UMAP and DM for each sample from AUC values
## Diffusion map and pseudotime based on UMAP coords of GRN

suppressPackageStartupMessages({
  library(destiny)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(scater)
  library(AUCell)
  
})

range_y <- function(x){
  mx <- quantile(x,0.99)
  mn <- quantile(x,.1)
  x[x>mx] <- mx
  x[x<mn] <- mn
  a <- (mx-mn)/(max(x)-min(x))
  b <- mx-(a*max(x))
  d <- round((a*x)+b, digits = 3)
  return(d)
  
}
range_x <- function(x){
  mx <- quantile(x,0.99)
  mn <- quantile(x,.1)
  x[x>mx] <- mx
  x[x<mn] <- mn
  return(x)
  
}
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
###DIRs
DIR_GRNind <- "SCENIC/GRNind/"
DIR_PS <- "~/MBsnANA/ATACana/ATAC_ArchR/Piyush_share/"
DIR_GRNana <- "SCENIC/GRNana/"
####
Sample <- commandSample(trailingOnly = TRUE)
msn <- gsub("MSM","MSN",Sample)
#
# #############UMAP and diffusion map based on AUC values-##############
# ###load data
## integrated tumor plot data after cluster annotation
# load(paste0(DIR_GRNana,"plotdata_G34MBaucKNN.RData"))
# plot.data <- plot.data.com[plot.data.com$Batch==msn,]
# rm(plot.data.com)
## load TF-GRN AUC values
# load(paste0(DIR_GRNana,"pdG34MB_G34GRN_AUC.RData"))
# scr <- t(assay(ms))
# table(row.names(plot.data) %in% row.names(scr))
# scr <- scr[row.names(plot.data),]
# 
## UMAP
# set.seed(1)
# um <- uwot::umap(scr, min_dist = 0.5, metric = "cosine")
# 
# plot.data$aucUMAP1 <- um[,1]
# plot.data$aucUMAP2 <- um[,2]
# 
# values <- read.delim("Files/anncol.txt",row.names = 1)
# d <- sort(unique(plot.data$Annotation))
# valuesx <- values[d,]
# 
# label.d = plot.data %>% group_by(Annotation) %>% 
#   select(aucUMAP1, aucUMAP2) %>% summarize_all(median)
# #UMAP ANN w labels
# p <- ggplot(plot.data, aes(x=aucUMAP1, y=aucUMAP2 ,color=Annotation))+
#   scale_color_manual(values = valuesx)+
#   geom_point()+
#   geom_label_repel(aes(label = Annotation),size = 2.5, data = label.d, show.legend = FALSE)+
#   theme_classic()+
#   theme(axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title=element_blank(),
#         legend.position = "None")
# tiff(paste0(DIR_GRNind,"figures/",msn,"_UMAP_AUC.tiff"), width = 6, height = 5, units = "in", res = 300)
# print(p)
# dev.off()
# #UMAP ANN
# p <- ggplot(plot.data, aes(x=aucUMAP1, y=aucUMAP2 ,color=Annotation))+
#   scale_color_manual(values = valuesx)+
#   geom_point(size=0.5)+
#   theme_classic()+
#   theme(axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title=element_blank(),
#         legend.position = "None")
# tiff(paste0(DIR_GRNind,"figures/",msn,"_UMAP_AUC_ANN.tiff"), width = 4, height = 3, units = "in", res = 300)
# print(p)
# dev.off()
# values <- read.delim("Files/axiscol.txt",row.names = 1)
# d <- sort(unique(plot.data$Axis))
# valuesx <- values[d,]
# #UMAP ANN
# p <- ggplot(plot.data, aes(x=aucUMAP1, y=aucUMAP2 ,color=Axis))+
#   scale_color_manual(values = valuesx)+
#   geom_point(size=0.5)+
#   theme_classic()+
#   theme(axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title=element_blank(),
#         legend.position = "None")
# tiff(paste0(DIR_GRNind,"figures/",msn,"_UMAP_AUC_Axis.tiff"), width = 4, height = 3, units = "in", res = 300)
# print(p)
# dev.off()
# 
# save(plot.data, file = paste0(DIR_GRNind,msn,"_AUC_UMAP.RData"))
# 
# ##diffusion map
# dm <- DiffusionMap(um, n_pcs=NA)
# a <- data.frame(DC1 = eigenvectors(dm)[, 1], DC2 = eigenvectors(dm)[, 2])
# pd <- cbind(plot.data,a)
# ##based on cell cycle
# k1 <- order(pd$CC.score, decreasing = TRUE)
# tps <- k1[1:3]
# dpt <- DPT(dm, tips = tps)
# #get average dpt value from tips
# dptval <- dpt[tps,]
# dptval <- colMeans(dptval)
# pd$DPTval <- dptval
# #
# ##UMAP DPTval--------------
# p <- ggplot(pd, aes(x=aucUMAP1, y=aucUMAP2, color=DPTval)) +
#   geom_point(size=0.5)+
#   scale_color_viridis_c(option = "C")+
#   theme_void()+
#   theme(legend.position='none',
#         axis.title = element_blank(),
#         axis.text.y = element_blank(),
#         axis.line = element_blank(),
#         axis.ticks = element_blank())
# tiff( paste0(DIR_GRNind,"figures/",msn,"_UMAP_DPT.tiff"), width = 4, height = 3, units = "in", res = 300)
# print(p)
# dev.off()
# 
# 
# p <- ggplot(pd, aes(x=DC1, y=DC2, color=DPTval)) +
#   geom_point(size=0.5)+
#   scale_color_viridis_c(option = "C")+
#   theme_void()+
#   theme(legend.position='none',
#         axis.title = element_blank(),
#         axis.text.y = element_blank(),
#         axis.line = element_blank(),
#         axis.ticks = element_blank())
# 
# tiff( paste0(DIR_GRNind,"figures/",msn,"_DMAP_DPT.tiff"),  width = 4, height = 3, units = "in", res = 300)
# print(p)
# dev.off()
# values <- read.delim(paste0(DIR_GRNana,"anncol.txt"),row.names = 1)
# d <- sort(unique(pd$Annotation))
# valuesx <- values[d,]
# p <- ggplot(pd, aes(x=DC1, y=DC2, color=Annotation)) +
#   scale_color_manual(values = valuesx)+
#   geom_point(size=0.5)+
#   theme_void()+
#   theme(axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         legend.position = "None")
# 
# tiff( paste0(DIR_GRNind,"figures/",msn,"_DMAP_ANN.tiff"),  width = 4, height = 3, units = "in", res = 300)
# print(p)
# dev.off()
# values <- read.delim(paste0(DIR_GRNana,"axiscol.txt"),row.names = 1)
# d <- sort(unique(pd$Axis))
# valuesx <- values[d,]
# p <- ggplot(pd, aes(x=DC1, y=DC2, color=Axis)) +
#   scale_color_manual(values = valuesx)+
#   geom_point(size=0.5)+
#   theme_void()+
#   theme(axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         legend.position = "None")
# 
# tiff( paste0(DIR_GRNind,"figures/",msn,"_DMAP_AX.tiff"),  width = 4, height = 3, units = "in", res = 300)
# print(p)
# dev.off()
# # 
# ##remove imputed cells
# pdx <- pd[pd$Imputed=="No",]
# pdx$CC.score <- range_01(pdx$CC.score)
# #UMAP CC
# p <- ggplot(pdx %>% arrange(CC.score), aes(x=aucUMAP1, y=aucUMAP2 ,color=CC.score))+
#   scale_color_gradient(low = "white",high="red4")+
#   geom_point(size=0.5)+
#   theme_classic()+
#   theme(axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title=element_blank(),
#         legend.position = "None")
# tiff( paste0(DIR_GRNind,"figures/",msn,"_UMAP_CC.tiff"),  width = 4, height = 3, units = "in", res = 300)
# print(p)
# dev.off()
# #DC CC
# p <- ggplot(pdx %>% arrange(CC.score), aes(x=DC1, y=DC2 ,color=CC.score))+
#   scale_color_gradient(low = "white",high="red4")+
#   geom_point(size=0.5)+
#   theme_classic()+
#   theme(axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title=element_blank(),
#         legend.position = "None")
# tiff( paste0(DIR_GRNind,"figures/",msn,"_DC_CC.tiff"),  width = 4, height = 3, units = "in", res = 300)
# print(p)
# dev.off()
# 
# save(pd, file =paste0(DIR_GRNind,msn,"_AUC_DM_UMAP.RData"))
# #q()
# 
####AUC score on UMAP#######################
load(paste0(DIR_GRNind,msn,"_AUC_DM_UMAP.RData"))

# plotgrn <- function(x){
#   scrx <- scr[row.names(pd),]
#   Gene=scale(range_x(scrx[,x]))
#   Gene[Gene>2] <- 2
#   Gene[Gene<(-2)] <- (-2)
#   pdy <- cbind(pd,Gene=Gene )
#   p <- ggplot(pdy %>% arrange(Gene), aes(x=aucUMAP1, y=aucUMAP2 ,color=Gene))+
#     scale_color_gradientn(colors = c("deepskyblue4","beige","magenta3"),    
#                           breaks=c(-2,0,2),
#                           limits=c(-2, 2))+
#     geom_point(size=0.5)+
#     theme_classic()+
#     #labs(title=x)+
#     theme(axis.line = element_blank(),
#           axis.ticks = element_blank(),
#           axis.text = element_blank(),
#           axis.title=element_blank(),
#           legend.position = "None")
#   return(p)
# }
# genes <- c("CRX","TBR1","LHX1")
# for(x in genes){
#   tiff( paste0(DIR_GRNind,"figures/",msn,"_G34aucGRN_",x,".tiff"),  width = 4, height = 3, units = "in", res = 300)
#   print(plotgrn(x))
#   dev.off()
# }
# rm(x)
# x <- "TBR1"
# plotgrndc <- function(x){
#   scrx <- scr[row.names(pd),]
#   Gene=scale(range_x(scrx[,x]))
#   Gene[Gene>2] <- 2
#   Gene[Gene<(-2)] <- (-2)
#   pdy <- cbind(pd,Gene=Gene )
#   p <- ggplot(pdy %>% arrange(Gene), aes(x=DC1, y=DC2 ,color=Gene))+
#     scale_color_gradientn(colors = c("deepskyblue4","beige","magenta3"), 
#                           breaks=c(-2,0,2),
#                           limits=c(-2, 2))+
#     geom_point(size=0.5)+
#     theme_classic()+
#     #labs(title=x)+
#     theme(axis.line = element_blank(),
#           axis.ticks = element_blank(),
#           axis.text = element_blank(),
#           axis.title=element_blank(),
#           legend.position = "None")
#   return(p)
# }
# x <- "RAX2"
# tiff( paste0(DIR_GRN,"figures/",msn,"_G34aucGRN_",x,"_DC.tiff"),  width = 4, height = 3, units = "in", res = 300)
# plotgrndc(x)
# dev.off()
####gene score on UMAp#######################
DIR_RNA <- "RNA/"
load(paste0(DIR_RNA,msn,"scepro.RData"))
mdt <- logcounts(get(msn))
# 
# plotgenedc <- function(x){
#   pdy <- pd[!(is.na(pd$ind_cluster)),]
#   mdtx <- mdt[,row.names(pdy)]
#   Gene=scale(range_x(mdtx[x,]))
#   Gene[Gene>2] <- 2
#   Gene[Gene<(-2)] <- (-2)
#   pdy$Gene <- Gene
#   p <- ggplot(pdy %>% arrange(Gene), aes(x=DC1, y=DC2 ,color=Gene))+
#     scale_color_gradientn(colors = c("skyblue","grey93","orangered4"), 
#                           breaks=c(-2,0,2),
#                           limits=c(-2, 2))+
#     geom_point(size=0.5)+
#     theme_classic()+
#     labs(title = x)+
#     theme(axis.line = element_blank(),
#           axis.ticks = element_blank(),
#           axis.text = element_blank(),
#           axis.title=element_blank(),
#           legend.position = "None")
#   return(p)
# }
# y<- "MYC"
# tiff( paste0(DIR_GRN,"figures/",msn,"_DC_",y,".tiff"),  width = 4, height = 3, units = "in", res = 300)
# plotgenedc(y)
# dev.off()

plotgene <- function(x){
  pdy <- pd[!(is.na(pd$ind_cluster)),]
  mdtx <- mdt[,row.names(pdy)]
  Gene=scale(range_x(mdtx[x,]))
  Gene[Gene>2] <- 2
  Gene[Gene<(-2)] <- (-2)
  pdy$Gene <- Gene
  p <- ggplot(pdy %>% arrange(Gene), aes(x=aucUMAP1, y=aucUMAP2 ,color=Gene))+
    scale_color_gradientn(colors = c("skyblue","grey93","orangered4"), 
                          breaks=c(-2,0,2),
                          limits=c(-2, 2))+
    geom_point(size=0.5)+
    theme_classic()+
    #labs(title = x)+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title=element_blank(),
          legend.position = "None")
  return(p)
}
genes <- c("CRX","PAX6","EOMES")
# for(x in genes){
#   tiff( paste0(DIR_GRNind,"figures/",msn,"_UMAP_",x,".tiff"),  width = 4, height = 4, units = "in", res = 300)
#   print(plotgene(x))
#   dev.off()
# }
#rm(x)
pdy <- pd[!(is.na(pd$ind_cluster)),]
mdtx <- as.matrix(mdt[genes,row.names(pdy)])
corx <- cor(t(mdtx),t(mdtx))
myColor <- colorRampPalette(c("blue","white","red4"))(100)
myBreaks <- c(seq(-1, 0, length.out=ceiling(100/2) + 1), 
              seq(1/100, 1, length.out=floor(100/2)))
hp <- pheatmap::pheatmap(corx,
                         color=myColor, breaks=myBreaks,
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         border_color = NA, 
                         fontsize = 7)

tiff(paste0(DIR_GRNind,"figures/Corr_",msn,".tiff"), width = 6, height = 6, units = "in", res = 300)
print(hp)
dev.off()

####G34MB AUC score on boxplot#######################
load(paste0(DIR_GRNana,"pdG34MB_G34MBNMF_AUC.RData"))
scr <- t(assay(ms))
rm(plot.data.com,ms)

keep <- pd$Axis %in% c("PR","Precursors","UBC")
pdx <- pd[keep,]
scrx <- scr[row.names(pdx),]
pdx$NMF1 <- scrx[,1]
pdx$NMF2 <- scrx[,2]
pdx$NMF1 <- range_01(pdx$NMF1)
pdx$NMF2 <- range_01(pdx$NMF2)

pdx2 <- reshape2::melt(pdx[,c("Axis","NMF1","NMF2")])
pdx2$variable <- factor(pdx2$variable, levels =c("NMF2","NMF1") )
p <- ggplot(pdx2, aes(x=Axis,y=value, fill=variable))+
  geom_boxplot()+
  scale_fill_manual(values = c("goldenrod1","green4"))+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black',linewidth=1),
        axis.ticks = element_line(colour = 'black',linewidth=1),
        axis.text.y = element_text(face="bold",colour = "black",size=15),
        #axis.text.x = element_text(face="bold",colour = "black",size=8, angle = 90,vjust=0.5, hjust=1),
        axis.text.x = element_blank(),
        axis.title=element_blank(),
        legend.position = "None")+
  scale_y_continuous(breaks =c(0, 0.5, 1))
tiff( paste0(DIR_GRNind,"figures/",msn,"_G34MB_NMFsig.tiff"),  width = 4.2, height = 4, units = "in", res = 300)
print(p)
dev.off()
q()
