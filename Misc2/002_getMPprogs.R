## This script uses NMF factors obtained across samples and rank combinations.
## For each factor, top 100 genes ranked (decrearing order) by contribution to 
## that factor are assigned as a representative meta-gene signature of that factor

## Then Jaccard and Overlap similaity is calculated across all sample X ranks
## meta-gene sets. Overlap similairty matrix is used to obtain robust meta-gene
## gene-sets, which are then merged into meta-gene programs 

suppressPackageStartupMessages({
  library(Matrix)
})

## function for jaccard similarity and overlap similarity
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (round(intersection/union,4))
}
overlapS <- function(a, b) {
  intersection = length(intersect(a, b))
  minlen = min(length(a),length(b))
  return (round(intersection/minlen,4))
}

DIR <- "RNA/NMFana/"
setwd(DIR)
########## Obtain Similairty matrix across all factors ##########
## 1. use 'for' loop to load nmfs output for each sample and rank combination-----

genelist <- list()
geneimplist <- list() ##need this for importance score later for merging

## two columns, Sample id, NMF ranks
sams <- readLines("Files/MBlibs.txt")
nmf_ranks <- seq(3,18,by=3)
all_runs <- row.names(sams_rnks)<-paste0(sams_rnks$V1,"_",sams_rnks$V2)

##loop on sample name

for(sam in sams){
  ##loop on nmf ranks
  for(sel_rank in nmf_ranks){
    #load nmf data
    load(paste0("MB_",sam,"_NMFr",sel_rank,".RData"))
    nmfs <- colnames(rd_nmf)
    for(nmfid in 1:length(nmfs)){
      del <- data.frame(Score=sort(H[nmfid,], decreasing = TRUE)[1:100]) ## top 100
      del$Gene <- row.names(del)
      del$Rank <- seq(dim(del)[1])
      genelist[[paste0(sam,"_",sel_rank,":",nmfid)]]  <- row.names(del)
      geneimplist[[paste0(sam,"_",sel_rank,":",nmfid)]] <- del
      rm(del)
    }
    rm(nmfid,nmfs,rd_nmf,H,pd)
    
  }
  rm(sel_rank)
  
}
rm(sam)


## 2. similarity of NMF gene-sets across all NMF gene-sets----------
all_nmfs <- names(genelist)

## Overlap similarity, vectorized
OS_AvB <- outer(all_nmfs, all_nmfs, FUN = Vectorize(function(i, j) {
  overlapS(genelist[[i]], genelist[[j]])
}))
rownames(OS_AvB) <- colnames(OS_AvB) <- all_nmfs


## Jaccard similarity, vectorized
JS_AvB <- outer(all_nmfs, all_nmfs, FUN = Vectorize(function(i, j) {
  jaccard(genelist[[i]], genelist[[j]])
}))
rownames(JS_AvB) <- colnames(JS_AvB) <- all_nmfs

save(JS_AvB, geneimplist, genelist, file = "JS_AvBforNMFngene100.RData")
save(OS_AvB, geneimplist, genelist, file = "OS_AvBforNMFngene100.RData")
q()

##### Obtain robust gene-sets and then gene-programs #######
## 3. Load Similarity Matrix -----------
load("OS_AvBforNMFngene100.RData")

## create a meta-data data-frame
sams <- readLines("Files/MBlibs.txt")
nmf_ranks <- seq(3,18,by=3)
all_runs <- row.names(sams_rnks)<-paste0(sams_rnks$V1,"_",sams_rnks$V2)

nmf_ranks <- seq(3,18,by=3)
ann <- data.frame(nmfs=names(genelist),
                  Dataset=rep(sams,each=63),
                  Rank=rep(nmf_ranks,nmf_ranks))
row.names(ann) <- ann$nmfs
ann$nmfs <- NULL
head(ann)


to_keep <- c() ## this will store representative robust gene-sets across samples

for(sam in sams){
  
  ## First the aim is to obtain robust nmf gene-set per rank per sample across ranks
  ## Robustness is defined by overlap across other factors
  ## Per rank robust gene-set have greater overlap with gene-set across other ranks
  ## Than redundancy is removed within a rank and than redundancy across other 
  ## ranks for the same sample
  
  ## obtain nmfs associated with a sample
  all_nmfs <- row.names(ann)[ann$Dataset==sam]
  ## subset the OS matrix to sample
  OS_AvB_CT <- OS_AvB[all_nmfs,all_nmfs]
  dim(OS_AvB_CT)
  ## obtain all the ranks
  ranks <- paste0("_",nmf_ranks,":") 
  
  ## this is to obtain robust non-redundant nmfs for each rank per sample 
  to_keep_a <- c() ## list of robust gene-sets across ranks
  for(x in ranks){
    ## all nmfs for a rank 
    sel_nmfs <- grep(x,all_nmfs)
    ##extract OS of a rank with all the other ranks
    delOS <- OS_AvB_CT[sel_nmfs,-c(sel_nmfs)]
    ## max similarity for each ranks nmfs with other nmfs, sorted by max overlap
    delOS2 <- sort(apply(delOS, 1, max), decreasing = TRUE)
    ##only keep robust nmfs for rank, 
    ## these nmfs overlap with at least one another nmf from other rank, CO=20%
    delOS2 <- delOS2[delOS2>=0.2]
    selx <- names(delOS2) ##robust part 1
    
    ## Remove redundant nmfs within the same rank
    delOS3 <- OS_AvB_CT[selx,selx]
    delOS4 <- matrix(0, nrow = dim(delOS3)[1],ncol = dim(delOS3)[1])
    row.names(delOS4) <- colnames(delOS4) <- row.names(delOS3)
    ##create a matrix where overlaping nmfs are marked by 1, CO=33%
    delOS4[delOS3> 0.33] <- 1 
    ## here the diag should not be set to 0, as it will removed in code below
    
    ## Removing nmfs from a rank that are overlapped by selected
    ## robust nmf from the same rank, sorted by max overlap with other ranks
    sel_nmfsy <- c()
    while (length(selx)>0) {
      if(length(selx)>1){
        ## take the first one
        selz <- selx[1]
        ## add it to selected nmf
        sel_nmfsy <- c(sel_nmfsy,selz)
        ##extract OS for other nmfs which overlap with selected nmf
        delOS5 <- colnames(delOS4)[delOS4[selz,]==1]
        ##remove other nmfs from that overlap the selected nmf from the selx
        delOS4 <- delOS4[-c(which(row.names(delOS4) %in% delOS5)),-c(which(row.names(delOS4) %in% delOS5))]
        selx <- selx[-c(which(selx %in% delOS5))]
      }else {
        sel_nmfsy <- c(sel_nmfsy,selx)
        selx <- c()
      } 
    }
    
    to_keep_a <- c(to_keep_a,sel_nmfsy)
    rm(sel_nmfs,delOS,delOS2,selx,delOS3,delOS4,sel_nmfsy,delOS5)
  }
  rm(x)
  
  ## Removing redundancy across ranks of the same sample
  ## Factors are ordered by max overlap across other samples 
  ## This also partially removes sample specific robust factors
  ## OS for other samples
  keep2 <-row.names(ann)[ann$Dataset!=sam]
  OS_AvB_CTx <- OS_AvB[to_keep_a,keep2]
  
  row_max <- apply(OS_AvB_CTx,1,max)
  row_max <- sort(row_max,decreasing = TRUE) ## these are also sorted in ascending rank
  sely <- names(row_max)
  
  ##repeating above method
  delOS3 <- OS_AvB_CT[sely,sely]
  delOS4 <- matrix(0, nrow = dim(delOS3)[1],ncol = dim(delOS3)[1])
  row.names(delOS4) <- colnames(delOS4) <- row.names(delOS3)
  ##create a matrix where overlaping nmfs are marked by 1, CO= 50%
  delOS4[delOS3>=0.5] <- 1 
  
  ## now removing  nmfs from a sample that are overlapped by selected
  ## robust nmf ranked by overlap across other ranks
  sel_nmfsy <- c()
  while (length(sely)>0) {
    if(length(sely)>1){
      ## take the first one
      selz <- sely[1]
      ## add it to selected nmf
      sel_nmfsy <- c(sel_nmfsy,selz)
      ##extract OS for other nmfs which overlap with selected nmf
      delOS5 <- colnames(delOS4)[delOS4[selz,]==1]
      ##remove other nmfs from that overlap the selected nmf from the selx
      delOS4 <- delOS4[-c(which(row.names(delOS4) %in% delOS5)),-c(which(row.names(delOS4) %in% delOS5))]
      sely <- sely[-c(which(sely %in% delOS5))]
      #print("reminaing nmfs :")
      #print(selx)
    }else {
      sel_nmfsy <- c(sel_nmfsy,sely)
      sely <- c()
    } 
  }
  
  to_keep <- c(to_keep,sel_nmfsy)
  rm(all_nmfs, OS_AvB_CT, OS_AvB_CTx,ranks,to_keep_a)
  rm(row_max,sely,delOS3,delOS4,sel_nmfsy,delOS5)
}
rm(sam)

length(to_keep) ## these are filteres robust gene-sets across samples and ranks

## 4. Using Jaccard similairty to hierarchially cluster meta-gene sets into programs----------- 
library(pheatmap)
load("JS_AvBforNMFngene100.RData")
JS <- JS_AvB[to_keep,to_keep]

## obtain broad clusters
hp2 <-pheatmap(OS,clustering_method = 'ward.D2',
               show_rownames = FALSE,show_colnames = FALSE)

del2 <- hp2$tree_col
del3 <- del2$labels
del3 <- del3[del2$order]
cluster = data.frame(X=cutree(hp2$tree_row, k = 40)) # starting with 40 MPs
head(cluster)

ann$X <- "X"
annx <- ann[to_keep,]
annx$MP <- paste0("MP",cluster$X)
annx$X <- NULL
annx$Rank <- NULL
## save first broad clusters
write.table(annx, file = "ann_40wd2_nmfMP100g_intial.txt", sep = "\t", quote = FALSE)

## 5. Alliteratively merging MPs--------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
})
load("OS_AvBforNMFngene100.RData")
mp_genes <- list()

## intial
annx <- read.delim("ann_40wd2_nmfMP100g_intial.txt",row.names = 1)
annx$NewMP <- annx$MP
nmfmpclust <- unique(annx$MP)

# subsequent iterations
#annx <- read.delim("nmf100g_mpfix.txt",row.names = 1) ## iterations
#nmfmpclust <- unique(annx$NewMP)

for(mpx in nmfmpclust){
  #identify nmfs clustered together
  sel_nmfs <- row.names(annx)[annx$NewMP==mpx]
  ##subset the list
  geneimplist_sub <- list()
  for(nmfx in sel_nmfs){
    geneimplist_sub[[nmfx]] <- geneimplist[[nmfx]]
    
  }
  rm(nmfx)
  ##merge the obtaind list of data frames
  merged_gene_df <- bind_rows(geneimplist_sub) %>%
    group_by(Gene) %>%
    summarise(Rank = mean(Rank, na.rm = TRUE), 
              count = n(), .groups = "drop")
  ##first rank by mean importance of genes and then by number of time gene appears
  ## this will result in most frequent gene and between gene that has same frequency 
  ##the one with higher mean importance/ or rank
  merged_gene_df_ranked <- merged_gene_df %>%
    arrange(desc(count), Rank)
  mp_genes[[mpx]] <- merged_gene_df_ranked$Gene[1:100] ## new representative gene-set
  
  rm(merged_gene_df,merged_gene_df_ranked,geneimplist_sub,sel_nmfs)
}
rm(mpx)
length(mp_genes)
lengths(mp_genes)


## 6. checking merged MPs for overlap ---------

all_mps <- names(mp_genes)
OS2<- outer(all_mps, all_mps, FUN = Vectorize(function(i, j) {
  overlapS(mp_genes[[i]], mp_genes[[j]])
}))
rownames(OS2) <- colnames(OS2) <- all_mps

OS3 <- matrix(0,length(all_mps),length(all_mps))
row.names(OS3) <- colnames(OS3) <- all_mps
OS3[OS2>0.5] <- 1 ## cut-off
## heatmap
mb2 <- ((seq(1:101)-1)/100)
mc3 <- colorRampPalette(c("white","grey87","black"))(100)

hp2 <-pheatmap(OS3,
               #clustering_method = 'ward.D2',
               cluster_rows = FALSE,cluster_cols = FALSE,
               breaks = mb2, color = mc3, fontsize = 8)
pdf("NMF_OS_50.pdf", width = 8, height = 8, pointsize = 10)
print(hp2)
dev.off()

##save this list of TopicMP CREs
save(mp_genes, file = "NMFMPgene_26MPfix_100g.RData")
q()



