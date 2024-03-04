##This script gets combined RNA-ATAC data from scJoint run
suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggsci)
  library(ggrepel)
})
##DIRs
DIR_RNA <- "/Retina/RNA/"
DIR_ATAC <- "Retina/ATAC/"
DIR_scJoint <- "Retina/ATAC/scJoint/output/"


load(paste0(DIR_ATAC,"com_RNA_ATAC_pd.RData"))

df <- read.delim(paste0(DIR_scJoint,"Retina_df.txt")) ##remove space after name

plot.data <- cbind(plot.data,df)
##checking how many predicted are same
keep <- plot.data$Platform=="RNA"
table(plot.data[!keep,"Annotationlv2"]==plot.data[!keep,"predicted"]) ##all RNA ones are same
save(plot.data, file = paste0(DIR_ATAC,"com_RNA_ATAC_pd_scJ.RData"))
d <- length(unique(plot.data$predicted))
label.d = plot.data %>% group_by(predicted) %>% 
  select(UMAP1, UMAP2) %>% summarize_all(median)
#for plot
p <- ggplot(plot.data, aes(x=UMAP1, y=UMAP2, color=predicted))+
  geom_point(size=0.5)+
  scale_color_manual(values =colorRampPalette(pal_futurama()(12))(d))+
  #scale_color_manual(values =valuesx)+
  theme_classic()+
  geom_label_repel(aes(label = predicted),size = 2.5, data = label.d, show.legend = FALSE)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "None")+
  facet_wrap(~Platform)

tiff("Retinac.tiff", unit = "in",  width=12, height=6, res = 300)
print(p)
dev.off()
