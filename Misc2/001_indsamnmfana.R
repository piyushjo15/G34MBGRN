## This script obtains NMF factors for each tumor sample for multiple ranks

library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/path/to/python")
use_python("/path/to/python", required = TRUE)


suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(Matrix)
})

DIR_RNA <- "RNA/"

args <- commandArgs(trailingOnly = TRUE) 
sam <- args[1] ## sample id
rnk <- args[2] ## rank from 3 to 18 in increment of 3
rnk <- as.integer(as.numeric(rnk))
print(paste0("Performing NMF factoriziation for sample: ", sam, " and rank: ",rnk))

## scater/scran processed data for each sample
sce <- get(load(paste0(DIR_RNA,sam,"scepro.RData"))) 

## integrated tumor cell metadata
load("/path/to/Files/plotdata_G34MBaucKNN.RData")

## subset to a sample
pd <- plot.data.com[plot.data.com$Batch==sam,]

## combined tumor HVG, no ribosomal, mitochondiral or sex chromosome genes
load("topHVGG34new10xTUMnoSX_1500.RData") 
com.hvg <- readLines("/path/to/Files/topHVGG34new10xTUMnoSX_1500.txt") 


mdt <- t(logcounts(sce)[com.hvg,row.names(pd)]) ## transpose for NMF input

repl_python()
import sklearn.decomposition as sk
import numpy
import random

rnk=r.rnk
model = sk.NMF(n_components=rnk, init='nndsvda', max_iter=10000, random_state=0)
random.seed(456)
w = model.fit_transform(r.mdt)
h = model.components_

exit

rd_nmf= py$w
H = py$h

row.names(rd_nmf) <- row.names(pd)
colnames(rd_nmf) <- row.names(H) <- paste0(sam,"_",rnk,":",seq(dim(rd_nmf)[2]))
colnames(H) <- com.hvg

##save NMF factor, gene contribution to factors, sample metadata and HVG gene list
save(rd_nmf, pd, H, com.hvg,file = paste0("NMFana/MB_",sam,"_NMFr",rnk,".RData")) 

q()
