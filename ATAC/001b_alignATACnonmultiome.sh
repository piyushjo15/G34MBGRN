#!/bin/bash


##this script is for aligning ATAC data for non-multiome ATAC data
export PATH=$PATH:$HOME/path/to/CRatac # location of cellranger-atac

cellranger-atac count --id=${SCE} \
     # location of HG38 index for cellrangerATAC/ARC  
		--reference=$HOME/HB/index/CR_ARC_HG38 \
		--fastqs=$HOME/path/to/atac/fastqs/${SCE} \
		--sample=${SCE} 
