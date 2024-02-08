library(reticulate)
##python
Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")
#argument
suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(Matrix)
})

args = commandArgs(trailingOnly=TRUE) #sample id
load(paste0(args[1],"_postSoupXDX.RData")) #load output of diem, SoupX and decontX processing
pre_doublets <- colnames(sce) ##get raw counts

#scrublet----
mdt <- counts(sce)
mdt <- t(mdt)
# python code
repl_python()
import scrublet as scr
scrub = scr.Scrublet(r.mdt)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
##for MSN105 and MSN058
#predicted_doublets = scrub.call_doublets(threshold=0.8)
exit

sce$SCRscore <- py$doublet_scores
keep <- py$predicted_doublets
SCRcall <- rep("Singlet", length(keep))
SCRcall[keep] <- "Doublet"
sce$SCRcall<- SCRcall
sce <- subset(sce, ,SCRcall=="Singlet")
post_doublets <- colnames(sce)
doublets <- pre_doublets[!(pre_doublets %in%post_doublets )]
doublets <- gsub(paste0(args[1],"_"),"",doublets)
write(doublets,paste0(args[1],"_doublets.txt")) ##write identified doublets

save(sce, mdfr2,mdt.pd, file = paste0(args,"_postSCR.RData"))
q()
