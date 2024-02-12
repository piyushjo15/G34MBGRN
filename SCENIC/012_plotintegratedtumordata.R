suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggrepel)
  library(ggsci)
  library(plotly)
  
})

DIR_GRN <- "SCENIC/GRNana/"
load(paste0(DIR_GRN,"diff_G4MB_GRN_AUC.RData"))
plot.data.com <- cbind(plot.data.com,a)
#Axis colors
values <- read.delim("Files/axiscol.txt",row.names = 1)
d <- sort(unique(plot.data.com$Axis))
#cell-state colors
values <- read.delim("Files/anncol.txt",row.names = 1)
d <- sort(unique(plot.data.com$Annotation))
#subtype colors
values <- read.delim("Files/Subtype_G34.txt",row.names = 1)
d <- sort(unique(plot.data.com$Subtype))
valuesx <- values[d,]
#subgroup colors
valuesx <- c("goldenrod1","darkgreen")

##plotting diffusion map
plot_ly(data = plot.data.com, 
        x = ~DC1, y = ~DC2, z = ~DC3, 
        color = ~Annotation, 
        type = "scatter3d",
        colors = valuesx,
        #colors = c("grey","red4"),
        #colors = colorRampPalette(pal_simpsons()(16))(5),
        mode = "markers", 
        marker = list(size = 1, width=1), # controls size of points
        #text=~Axis, #This is that extra column we made earlier for which we will use for cell ID
        #hoverinfo="text",
        showlegend=FALSE) %>%
  layout(scene=list(
    xaxis =list(showgrid=FALSE, zeroline =FALSE),
    yaxis = list(showgrid=FALSE, zeroline =FALSE),
    zaxis = list(showgrid=FALSE, zeroline =FALSE),
    camera = list(eye = list(x = 1, y = 2, z = 3))))

##Cell cycle score
ggplot(plot.data.com %>% arrange(CC.score), aes(x=nfUMAP1, nfUMAP2, color=CC.score))+
  geom_point(size=0.3)+
  scale_color_gradient(high = "red",low = "white")+
  theme_classic()+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
##subtype
label.d = plot.data.com %>% group_by(Subtype) %>% 
  select(nfUMAP1, nfUMAP2) %>% summarize_all(median)

ggplot(plot.data.com, aes(x=nfUMAP1, nfUMAP2, color=Subtype))+
  geom_point(size=0.3)+
  theme_classic()+
  scale_color_manual(values = valuesx)+
  geom_label_repel(aes(label = Subtype),size = 2.5, data = label.d, show.legend = FALSE)+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
##Group
label.d = plot.data.com %>% group_by(Group) %>% 
  select(nfUMAP1, nfUMAP2) %>% summarize_all(median)

ggplot(plot.data.com, aes(x=nfUMAP1, nfUMAP2, color=Group))+
  geom_point(size=0.3)+
  theme_classic()+
  scale_color_manual(values = valuesx)+
  geom_label_repel(aes(label = Group),size = 2.5, data = label.d, show.legend = FALSE)+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

##clusters
d <- length(unique(plot.data.com$KNNcl))
label.d = plot.data.com %>% group_by(KNNcl) %>% 
  select(nfUMAP1, nfUMAP2) %>% summarize_all(median)

ggplot(plot.data.com, aes(x=nfUMAP1, nfUMAP2, color=KNNcl))+
  geom_point(size=0.3)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(pal_igv()(52))(d))+
  geom_label_repel(aes(label = KNNcl),size = 2.5, data = label.d, show.legend = FALSE)+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())