library(reticulate)

##python
Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")

##This script is to to obtain diffusion component for G34MB cells using NMF factors
## obtained from on AUC values for each cells using the combined GRNs
library(destiny)
library(uwot)
##diffusion ----
DIR_GRN <- "SCENIC/GRNana/"
#
#using NMF------
load(paste0(DIR_GRN,"pdG34MB_G34GRN_AUC_NMF.RData"))

load(paste0(DIR_GRN,"plotdata_G34MBaucKNN.RData"))
# diffusion
dm <- DiffusionMap(rd, k=25,  n_pcs = NA, 
                   distance = "cosine",
                   suppress_dpt=TRUE)
a <- data.frame(DC1 = eigenvectors(dm)[, 1], DC2 = eigenvectors(dm)[, 2],
                DC3 = eigenvectors(dm)[, 3])



save(plot.data.com, a,  file =paste0(DIR_GRN,"diff_G4MB_GRN_AUC.RData"))
q()