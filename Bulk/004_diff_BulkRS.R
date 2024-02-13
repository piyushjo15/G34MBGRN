##This script is to to obtain diffusion component for G34MB cells using NMF factors
## obtained from on AUC/MS values for each cells using the combined GRNs
library(destiny)
library(uwot)
##diffusion ----
load("NMF_G34bulkRS.RData")
# diffusion
set.seed(123)
dm <- DiffusionMap(rd_8, k=5,  n_pcs = NA,
                   suppress_dpt=TRUE, n_eigs=3)
a <- data.frame(DC1 = eigenvectors(dm)[, 1], DC2 = eigenvectors(dm)[, 2],
                DC3 = eigenvectors(dm)[, 3])


save(metax, a,  file= "diff_G4MB_NMF.RData")
q()
load("diff_G4MB_NMF.RData")
d <- sort(unique(metax$Subtype_M))
subtype <- read.delim("Subtype_G34.txt",row.names = 1)
subtypex <- subtype[d,1]

metaxd <- cbind(metax,a)
p <- ggplot(metaxd, aes(x=DC1, DC3, color=Subtype_M))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values = subtypex)+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
tiff( "G34_MB_diff.tiff", width = 6, height = 6, units = "in", res = 300)
print(p)
dev.off()
p <- ggplot(metaxd, aes(x=DC1, DC3, color=Subgroup_M))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values = groupx)+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
tiff( "G34_MB_diff_group.tiff", width = 6, height = 6, units = "in", res = 300)
print(p)
dev.off()
p <- ggplot(metaxd, aes(x=DC1, DC3, color=Batch))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values = dsc)+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
tiff( "G34_MB_diff_batch.tiff", width = 6, height = 6, units = "in", res = 300)
print(p)
dev.off()
load("pdMB_WGCNA_AUC.RData")
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

subs <- row.names(subtype)

rdx <- apply(rd_8, 2, range_01)
library(dplyr)
library(tidyr)
plotNMFDC <- function(x){
  pd <- cbind(metaxd[,c("DC1","DC3","Subgroup_M","Subtype_M")], NMF=rd_8[,x])
  p <- ggplot(pd %>% arrange(NMF), aes(x=DC1, y=DC3, color=NMF))+
    geom_point()+
    scale_color_gradient(low="white",high="red4")+
    theme_classic()+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title=element_blank(),
          legend.position = "None")
  return(p)
}
for(i in 1:dim(rd_8)[2]){
  tiff(paste0("Signautre_NMF_",i,"_DC.tiff"), width = 4, height = 3, units = "in", res = 300)
  print(plotNMFDC(i))
  dev.off()
}
sel2 <- c("green","lightyellow","turquoise","black",  "pink","brown")
sel2 <-colnames(scr)
load("pdMB_WGCNA_AUC.RData")
plotWGCNADC <- function(x){
  pd <- cbind(metaxd[,c("DC1","DC3","Subgroup_M","Subtype_M")], NMF=range_01(scr[row.names(metaxd),x]))
  p <- ggplot(pd %>% arrange(NMF), aes(x=DC1, y=DC3, color=NMF))+
    geom_point()+
    scale_color_gradient(low="white",high="red4")+
    theme_classic()+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title=element_blank(),
          legend.position = "None")
  return(p)
}
for(i in sel2){
  tiff(paste0("Signautre_WGCNA_",i,"_DC.tiff"), width = 4, height = 3, units = "in", res = 300)
  print(plotWGCNADC(i))
  dev.off()
}
plotWGCNADC(24)
