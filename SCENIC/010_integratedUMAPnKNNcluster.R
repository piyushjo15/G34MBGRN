##using module scores to perform factorization, get integrated UMAp representation
## and then call clusters
library(reticulate)
##python
Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")
#function for  return nearest neighbour less than radius and max neighbors
suppressPackageStartupMessages({
  library(Matrix)
  library(uwot)
  library(scater)
  library(AUCell)
})
DIR_GRN <- "SCENIC/GRNana/"

## loading the AUC score for the combined single-cell tumor data-------
load(paste0(DIR_GRN,"pdG34MB_G34GRN_AUC.RData") )

ms1 <- t(assay(ms)) ##convert into matrix and transpose
ms <- ms1
###python run for NMF factorization-----
repl_python()
import sklearn.decomposition as sk
import numpy as np
import random

model = sk.NMF(n_components=25,init="nndsvd", max_iter=10000, random_state=0, tol=1e-5)
random.seed(456)
model.fit(r.ms)
wn = model.transform(r.ms)

exit
rd <- py$wn
set.seed(134)
del1 <- umap(rd, n_components = 2, min_dist = 0.4, metric = "cosine", n_neighbors = 25)
colnames(del1) <-  c("nfUMAP1","nfUMAP2")
set.seed(134)
plot.data.com <- cbind(plot.data.com,del1)
save(rd, plot.data.com, file=paste0(DIR_GRN,"pdG34MB_G34GRN_AUC_NMF.RData"))

##python run for KNN-leiden clustering------

repl_python()
import numpy as np
import umap
import random
import igraph as ig
import leidenalg as la
import networkx as nx
import scipy.sparse as ss
import sklearn.neighbors as sk

wn=r.rd

##using NMF factors to identify nearest neighbors using sklearn neighbors graph
## This uses data to give a neighborhood graph
random.seed(455)
a = sk.kneighbors_graph(wn,n_neighbors=26, metric="cosine", n_jobs=3, include_self=True) 

## creating a networkx graph
nxg = nx.from_scipy_sparse_matrix(a) 
#converting to igraph
g_rnn = ig.Graph(len(nxg), list(zip(*list(zip(*nx.to_edgelist(nxg)))[:2])), directed=False)
del(nxg,a)

#getting edgelist
ela = g_rnn.get_edgelist()
## gettign weight list
wj = g_rnn.similarity_jaccard(pairs=ela) #for edge list pairs
#getting clustering
part = la.find_partition(g_rnn, la.ModularityVertexPartition, weights=wj, n_iterations= -1, max_comm_size=4000, seed=12)

# #another approach for resolution
# opt = la.Optimiser()
# opt.consider_empty_community = True
# opt.max_comm_size = 1000
# opt.set_rng_seed(123)
# part = la.RBConfigurationVertexPartition(g_rnn,  resolution_parameter = 0.3, weights=wj)
# opt.optimise_partition(part, n_iterations=-1)

k6 = part.modularity
#clusters
k7 = np.array(part.membership)

exit
k10 <- py$k7

print(py$k6)
#assigning clustering-
plot.data.com$KNN_cl <- paste0("msKNN_",k10)

save(plot.data.com, file = paste0(DIR_GRN,"plotdata_G34MBaucKNN.RData"))
q()
#-----
library(tidyr)
library(dplyr)
library(ggrepel)
library(ggsci)
load(paste0(DIR_GRN,"plotdata_G34MBaucKNN_0909.RData"))
sub <- read.delim(paste0(DIR_SCP,"Subtype_G34.txt"), row.names = 1)

d <- sort(unique(plot.data.com$Subtype))
valuesx <- sub[d,"Color"]

ggplot(plot.data.com, aes(DC2, DC3, color=Subtype))+
  geom_point(size=0.3, alpha=0.5)+
  scale_color_manual(values = valuesx)+
  theme_classic()+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
  facet_wrap(~Subtype,ncol=4)

ggplot(plot.data.com %>% arrange(CC.score), aes(x=nfUMAP1, nfUMAP2, color=CC.score))+
  geom_point(size=0.3)+
  scale_color_gradient(high = "red",low = "white")+
  theme_classic()+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())


d <- length(unique(plot.data.com$Subtype))
label.d = plot.data.com %>% group_by(Subtype) %>% 
  select(nfUMAP1b, nfUMAP2b) %>% summarize_all(median)

ggplot(plot.data.com, aes(x=nfUMAP1b, nfUMAP2b, color=Subtype))+
  geom_point(size=0.3)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(pal_simpsons()(14))(d))+
  geom_label_repel(aes(label = Subtype),size = 2.5, data = label.d, show.legend = FALSE)+
    theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

d <- length(unique(plot.data.com$Subtype))
label.d = plot.data.com %>% group_by(Subtype) %>% 
  select(nfUMAP1, nfUMAP2) %>% summarize_all(median)

ggplot(plot.data.com, aes(x=nfUMAP1, nfUMAP2, color=Subtype))+
  geom_point(size=0.3)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(pal_simpsons()(14))(d))+
  geom_label_repel(aes(label = Subtype),size = 2.5, data = label.d, show.legend = FALSE)+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

ggplot(plot.data.com, aes(x=nfUMAP1b, nfUMAP2b, color=ATAC))+
  geom_point(size=0.3)+
  theme_classic()+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

pd <- plot.data.com
load("pdG34MB_TFselClGRN_MSv2.RData")
sel <- read.delim("TFselClGRNs_sel.txt")
sel <- sel$GRNs
pdms <- plot.data.com[row.names(pd),sel]
ms2 <- apply(pdms, 2, zscore_c)
plotMS <- function(x){
  pdx <- cbind(pd,ms2[,x])
  colnames(pdx)[dim(pdx)[2]] <- "MS"
  p <- ggplot(pdx %>% arrange(MS), aes(x=nfUMAP1, nfUMAP2, color=MS))+
    geom_point(size=0.3)+
    scale_color_gradient2(high = "magenta3",mid = "beige",low = "deepskyblue4")+
    theme_classic()+
    theme(legend.position='none',
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank())
  return(p)
}
plotMS("MXD3")

knns <- data.frame(table(pd$KNN_cl))
keep <- knns$Freq>300
cl_lv <- as.character(knns$Var1[keep])

msd2 <- c()

for(x in cl_lv){
  keep1 <- pd$KNN_cl==x
  msx <- ms2[keep1,]
  del <- colMeans(msx)
  msd2 <- rbind(msd2,del)
  rm(keep1,msx,del)
}
rm(x)
rownames(msd2) <- cl_lv
msd2[1:4,1:4]
hp <- pheatmap::pheatmap(t(msd2[cl_lv,sel]), clustering_method = "ward.D2",
                         color=myColor, breaks=myBreaks,
                         cluster_rows = FALSE,  
                         fontsize = 8, 
                         border_color = NA)
myBreaks <- c(seq(-2, 0, length.out=ceiling(100/2) + 1), 
              seq(2/100,2, length.out=floor(100/2)))
myColor <- colorRampPalette(c("deepskyblue4","beige","magenta3"))(100)
