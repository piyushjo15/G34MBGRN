#getting metagene signature
suppressPackageStartupMessages({
  library(ggplot2)
  library(enrichR)
})
load("NMF_G34bulkRS.RData")

###G34MB meta signature
G34MB <- list()
G34MB[["NMF_1"]] <- tvg[order(h1[1,], decreasing = TRUE)[1:100]]
G34MB[["NMF_2"]] <- tvg[order(h1[2,], decreasing = TRUE)[1:100]]

###G34MB meta signature
G34MBSub <- list()
cl <- row.names(h2)
for(x in cl){
  G34MBSub[[x]] <- tvg[order(h2[x,], decreasing = TRUE)[1:100]]
}

terms <- c()
for(x in row.names(h1)){
  genes <- G34MB[[x]]
  ks <- enrichr(genes = genes,databases = "GO_Biological_Process_2023")
  ks <- ks$GO_Biological_Process_2023
  ks <- ks[order(ks$P.value),]
  ks <- ks[1:10,c("Term","Adjusted.P.value","Genes")]
  ks$Cluster <- x
  terms <- rbind(terms,ks)
  rm(genes,ks)
}
G34MB_BP <- terms
rm(x)

terms <- c()
for(x in cl){
  genes <- G34MBSub[[x]]
  ks <- enrichr(genes = genes,databases = "GO_Biological_Process_2023")
  ks <- ks$GO_Biological_Process_2023
  ks <- ks[order(ks$P.value),]
  ks <- ks[1:10,c("Term","Adjusted.P.value","Genes")]
  ks$Cluster <- x
  terms <- rbind(terms,ks)
  rm(genes,ks)
}
G34MBSub_BP <- terms
save(G34MB,G34MBSub, file = "G34MB_NMF_metagene.RData")
save(G34MB_BP,G34MBSub_BP, file = "G34MB_NMF_GO.RData")
q()
##jaccard index-----
library(RColorBrewer)
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
load("G34MB_NMF_metagene.RData")

sams <- names(G34MB)
G34MB_JS <- matrix(0, nrow = length(sams),ncol = length(sams))

for(i in 1:length(sams)){
  for (j in 1:length(sams)){
    G34MB_JS[i,j] <- jaccard(G34MB[[i]],G34MB[[j]])
  }
}
colnames(G34MB_JS) <- row.names(G34MB_JS) <- sams

mc <- colorRampPalette(c("white",brewer.pal(9,"BuPu")))(100)
mb <- ((seq(1:101)-1)/100)
hp <- pheatmap::pheatmap(G34MB_JS, 
                         cluster_rows = FALSE, cluster_col = FALSE,
                         breaks = mb, color = mc)
tiff("NMF_G34_JS.tiff", width = 3, height = 3.5, units = "in", res = 300)
print(hp)
dev.off()

sams <- names(G34MBSub)
G34MBs_JS <- matrix(0, nrow = length(sams),ncol = length(sams))

for(i in 1:length(sams)){
  for (j in 1:length(sams)){
    G34MBs_JS[i,j] <- jaccard(G34MBSub[[i]],G34MBSub[[j]])
  }
}
colnames(G34MBs_JS) <- row.names(G34MBs_JS) <- sams
sigsx2 <- paste0("NMF_",c(3,8,2,7,5,6,4,1))

mc <- colorRampPalette(c("white",brewer.pal(9,"BuPu")))(100)
mb <- ((seq(1:101)-1)/100)
hp <- pheatmap::pheatmap(G34MBs_JS[sigsx2,sigsx2], 
                         cluster_rows = FALSE, cluster_col = FALSE,
                         breaks = mb, color = mc)
tiff("NMF_G34sub_JS.tiff", width = 6, height = 6.5, units = "in", res = 300)
print(hp)
dev.off()