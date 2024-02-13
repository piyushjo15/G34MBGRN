## This script covnerts the GRN obtained from SCENIC+ analysis and coverts it 
## into python dictory for input to Celloracle
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "~/.conda/envs/pjexv4/bin/python3.7")
use_python("~/.conda/envs/pjexv4/bin/python3.7")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
##Dir
DIR_GRN <- "SCENIC/scenicplus/"
DIR_GRNind <- "SCENIC/RNA/"
Sample <- commandSample(trailingOnly = TRUE)

##making disctionary--------
msn <- gsub("MSM","MSN",Sample)
load(paste0(DIR_GRN,Sample,"/output/metadata/",msn,"_gmt.RData")) ##get identified GRNs component per sample
repl_python()
import os
import pickle

Sample=r.Sample
Sample2=r.msn
dataDir = os.path.join('/home/p541i/MBsnANA/ATACana/SCENIC/RNA/',Sample+'/')

your_dict=r.regs

# Save the dictionary to a file
with open(dataDir+'TFGRN.pkl', 'wb') as file:
  pickle.dump(your_dict, file)

exit

q()