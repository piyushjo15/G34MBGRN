## This script combines GRNs from different stages to obtain a conserved GRN 
## for Retina
sams <- read.delim("Retina_SCENIC_sams.txt", header = FALSE)
sams <- sams$V1
DIR_GRN <- "Retina/ATAC/SCENIC/scenicplus/"

com <- c()
for(x in sams){
  load(paste0(DIR_GRN,x,"/output/",x,"_top3GRN_AUC.RData")) 
  com <- c(com,all_reg)
  rm(all_reg)
  
}
rm(x)
df <- data.frame(table(com))
df <- df[df$Freq>1,]
regs_com <- as.character(df$com)

ADJ <- list()

for(x in sams){
  load(paste0(DIR_GRN,x,"/output/",x,"_gmt.RData"))
  regx <- names(regs)
  regx <- intersect(regx,regs_com)
  GS <- list()
  for(y in regx){
    GS[[y]] <- regs[[y]]
  }
  ADJ[[x]] <- GS
  
  rm(GS, regx,y, regs)
}
rm(x)

regs <- list()
for( y in regs_com){
  sel <- c()
  for(x in sams){
    del <- ADJ[[x]]
    tfc <- names(del)
    if(y %in% tfc) {
      del2 <- del[[y]]
      sel <- c(sel,del2)
      rm(del2)
    }
    rm(del, tfc)
  }
  rm(x)
  
  if(length(sel)>1){
    ##removing lowly represented genes
    tb <- data.frame(table(sel))
    co <- max(round(max(tb$Freq)/5)-1,2)
    keep3 <- tb$Freq>co
    tb <- tb[keep3,]
    regs[[y]] <- as.character(tb$sel)
    rm(keep3)
  }
  rm(sel)
}
rm(y)

for(y in regs_com){
  del <- length(regs[[y]])
  if(del<15){regs[[y]] <- NULL}
}
rm(y)
lengths(regs)
save(regs, file = paste0(DIR_GRN,"TFselcl_GRN_Retina.RData"))
q()
