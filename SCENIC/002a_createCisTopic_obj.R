## This script converts ATAC data from ArchR into cisTopic format
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")
library(ArchR)

# setwd
args <- commandArgs(TRUE)
Sample <- args[1]
RNA_Sample <- gsub("MSM","MSN",Sample)

#Next, I load the reference genome
addArchRGenome("hg38")

DIR_RNA <- "RNA/"
DIR_ATAC <- "ATAC/"

#Load the ATAC project
proj <- loadArchRProject(path = paste0(DIR_ATAC,Sample))


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


##########remove normal cell clusters!! ############
#load the cells to be removed

remove_cells <- readLines(paste0(DIR_ATAC,"normal_cells/",Sample,"_normal_cells.txt"))
count_df <- count_df[,colnames(count_df) %ni% remove_cells]

print("Number of normal cells")
length(remove_cells)

print("Number of cells after removing Normal cells")
dim(count_df)[2]


# Define the folder path and name
folder_path <- paste0("SCENIC/scATAC/",Sample,"/")
#get cell_metadata
metadata <- as.data.frame(getCellColData(proj))
metadata <- metadata[rownames(metadata) %ni% remove_cells,]

a <- rownames(metadata)
b <- paste(a,"___cisTopic",sep = "")
rownames(metadata) <- b

head(metadata)
#use python to generate cisTopic object ----------------
repl_python()

import os
import pandas as pd
import numpy
import warnings
import pycisTopic
from scipy import io

Sample = r.Sample
projDir = os.path.join("SCENIC/scATAC",Sample)
tmpDir = 'TMP/'

#load countmatrix from R
count = r.count_df
cell_data = r.metadata

#Create cisTopic object

from pycisTopic.cistopic_class import *

count_matrix=count
path_to_blacklist='Files/hg38-blacklist.v2.bed'

path_to_fragments= os.path.join("ATAC/",Sample,"outs/atac_fragments.tsv.gz")
cistopic_obj = create_cistopic_object( fragment_matrix=count_matrix,
                                      path_to_blacklist=path_to_blacklist,
                                      path_to_fragments= path_to_fragments)

cistopic_obj.add_cell_data(cell_data)
print(cistopic_obj)
#Save
import pickle
pickle.dump(cistopic_obj,
            open(os.path.join(projDir, 'cistopic_obj.pkl'), 'wb'))
            
#Save Backup
import pickle
pickle.dump(cistopic_obj,
            open(os.path.join(projDir, 'cistopic_obj_Backup.pkl'), 'wb'))


exit

print(paste0("Generated cistopic object for: ", Sample))
q()
