#!/bin/bash


##this script is for aligning ATAC data for 10x multiome experiments

export PATH=$PATH:$HOME/path/to/CR_ARC # location of cellranger-arc

cellranger-arc count --id=${SCE} \
      # location of HG38 index for cellrangerATAC/ARC 
		--reference=$HOME/HB/index/CR_ARC_HG38 \
		## csv file for location of RNA and ATAC fastqs for multiome experiment
		--libraries=$HOME/MBsnANA/scripts/snMULTI/${SCE}.csv 
		
