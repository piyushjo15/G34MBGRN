
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")
library(ArchR)
addArchRGenome("hg38")

#load RNA SCE
Sample <- commandArgs(trailingOnly = TRUE)

DIR_RNA <- "Retina/RNA/"
DIR_ATAC <- "Retina/ATAC/ArchR/"
DIR_OUT <- "Retina/SCENIC/"
setwd(DIR_OUT)


#Load the ATAC project
sample_path  <- paste0(DIR_ATAC,Sample)

proj <- loadArchRProject(path = sample_path)

#retrieve the PeakMatrix
Count <- getMatrixFromProject(proj, useMatrix = "PeakMatrix", binarize = F)
count_mat <- assay(Count)
class(count_mat)

count_df <- as.data.frame(as.matrix(count_mat))
class(count_df)

cell_names <- as.list(colnames(count_mat))
names(cell_names) <- colnames(count_mat)

#add region names to the rows of the dataframe
region_names <- as.data.frame(Count@rowRanges)
region_names_vector <- paste0(region_names$seqnames,":",region_names$start,"-",region_names$end)

rownames(count_df) <- region_names_vector


#get cell_metadata
load(paste0(DIR_ATAC,Sample,"/",Sample,"_metadata.RData"))
cls <- data.frame(table(plot.data_ATAC$predictedGroup))
keep <- cls$Freq<5
cls_rem <- as.character(cls$Var1[keep])
## remove annotated cells with less than 5 cells as these will cause 
## problem in scenicplus object creation
plot.data_ATAC <- plot.data_ATAC[!(plot.data_ATAC$predictedGroup %in% cls_rem),]
count_df <- count_df[,row.names(plot.data_ATAC)] ##subejct counts
#metadata <- as.data.frame(getCellColData(proj))
#metadata <- metadata[rownames(metadata) %ni% remove_cells,]

a <- rownames(plot.data_ATAC)
b <- paste(a,"___cisTopic",sep = "")
rownames(plot.data_ATAC) <- b

head(plot.data_ATAC)
metadata <- plot.data_ATAC ## dot cannot be read properly in python
# write.csv(metadata, 'Metadata.csv',
#           quote = FALSE,row.names = TRUE)

#use python to generate cisTopic object ----------------
repl_python()

import os
import pandas as pd
import numpy
import warnings
import pycisTopic
from scipy import io
import pickle

Sample = r.Sample
projDir = os.path.join("Retina/SCENIC/",Sample)
tmpDir = 'TMP/'

#load countmatrix from R
count = r.count_df
cell_data = r.metadata

#Create cisTopic object

from pycisTopic.cistopic_class import *

count_matrix=count
path_to_blacklist='Files/hg38-blacklist.v2_noncanchr.bed'

path_to_fragments= os.path.join(" Retina/ATAC/",Sample,"_fragments.tsv.gz")
cistopic_obj = create_cistopic_object( fragment_matrix=count_matrix,
                                      path_to_blacklist=path_to_blacklist,
                                      path_to_fragments= path_to_fragments)

cistopic_obj.add_cell_data(cell_data)
print(cistopic_obj)
#Save
pickle.dump(cistopic_obj,
            open(os.path.join(projDir, 'cistopic_obj.pkl'), 'wb'))
            
#Save Backup
pickle.dump(cistopic_obj,
            open(os.path.join(projDir, 'cistopic_obj_Backup.pkl'), 'wb'))


exit

print(paste0("Generated cistopic object for: ", Sample))
q()
