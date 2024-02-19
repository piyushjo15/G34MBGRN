## Script for UMAP, diffusionmap and slingshot for cerebellar UBC lineage
suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(uwot)
})


###using liger corrected GCUBC##########
load("GCUBCumapLiger5.RData")
#subsetting to UB Clineage
st <-c("progenitor_RL_early", "progenitor_RL","GCP/UBCP","UBC_diff","UBC_Trpc3","UBC_Hcrtr2") 
keep <- plot.data$subtype %in% st
pdUBC <- plot.data[keep,]
rd <- rd_als_cor[row.names(pdUBC),]
set.seed(134)
um <- umap(rd, min_dist = 0.8, metric = "cosine")
pdUBC$UMAP1 <- um[,1]
pdUBC$UMAP2 <- um[,2]
library(destiny)

dm <- DiffusionMap(rd, n_pcs=NA)
a <- data.frame(DC1 = eigenvectors(dm)[, 1], DC2 = eigenvectors(dm)[, 2])
pdUBC <- cbind(pdUBC,a)
##based on cell cycle
k1 <- order(pdUBC$CC.score, decreasing = TRUE)
tps <- k1[1:3]

dpt <- DPT(dm, tips = tps)
#get average dpt value from tips
dptval <- dpt[tps,]
dptval <- colMeans(dptval)
pdUBC$DPTval <- dptval

##plotting-----
st <-c("progenitor_RL_early", "progenitor_RL","GCP/UBCP","UBC_diff","UBC_Trpc3","UBC_Hcrtr2")
valuesx1 #color scale for subtype

dls <- c("progenitor","GCP/UBCP","UBC_diff","UBC_defined")
valuesx1 #color scale for dev_state
##
pdUBC$subtype <- factor(pdUBC$subtype, levels = st)
pdUBC$dev_state <- factor(pdUBC$dev_state, levels = dls)

p <- ggplot(data=pdUBC, aes(x=UMAP1, y=UMAP2, color=subtype)) +
  geom_point(size=.5 )+
  scale_color_manual(values=valuesx1)+
  theme_classic()+
  theme(legend.position = "None",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank())
tiff("UBC_UMAP_ST.tiff", width = 6, height = 4, units = "in", res = 300)
print(p)
dev.off()

p <- ggplot(data=pdUBC, aes(x=UMAP1, y=UMAP2, color=dev_state)) +
  geom_point(size=.5 )+
  scale_color_manual(values=valuesx2)+
  theme_classic()+
  theme(legend.position='none',
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank())
tiff("UBC_UMAP_DS.tiff", width = 6, height = 4, units = "in", res = 300)
print(p)
dev.off()
p <- ggplot(data=pdUBC, aes(x=UMAP1, y=UMAP2, color=DPTval)) +
  geom_point(size=.5 )+
  scale_color_viridis_c(option="B", direction = -1)+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
##slingshot-----
library(slingshot)
cl <- as.character(pdUBC$subtype)
u_na <- pdUBC[,c("UMAP1","UMAP2")]
sl <- slingshot(u_na, cl, start.clus = "progenitor_RL_early",
                end.clus="UBC_Hcrtr2", reducedDim=NULL)

library(ggridges)
pdUBC$Rank <- sl@metadata$curves$Lineage1$lambda
p1 <- ggplot(pdUBC, aes(x=Rank, y=subtype, fill=subtype) ) +
  geom_density_ridges()+
  scale_fill_manual(values=valuesx1)+
  theme_classic()+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
tiff("UBC_PST_ST.tiff", width = 6, height = 4, units = "in", res = 300)
print(p1)
dev.off()
p2 <- ggplot(pdUBC, aes(x=Rank, y=dev_state, fill=dev_state) ) +
  geom_density_ridges()+
  scale_fill_manual(values=valuesx2)+
  theme_classic()+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
tiff("UBC_PST_DS.tiff", width = 6, height = 4, units = "in", res = 300)
print(p2)
dev.off()
p3 <- ggplot(data=pdUBC, aes(x=UMAP1, y=UMAP2, color=Rank)) +
  geom_point(size=.5 )+
  scale_color_viridis_c(option = "B")+
  theme_classic()+
  theme(legend.position = "None",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
tiff("UBC_PST_UMAP.tiff", width = 6, height = 4, units = "in", res = 300)
print(p3)
dev.off()
save(pdUBC,sl, file = "UBC_slingshot_lig5.RData")
q()
