#!/bin/bash

##activate conda env for pyscenic run
module load conda
source activate pyscenic
echo "Processing sample s..../n"
echo $SCE
##using loom file generated from prepforSCENIC to get adjcancey matrcex
pyscenic grn ${SCE}_nn.loom hg_tfs.txt -m genie3 -o adj${SCE}.csv --num_workers 8 --seed 10
#add correlation
pyscenic add_cor adj${SCE}.csv ${SCE}_nn.loom \
    --output cor${SCE}.tsv \
    --mask_dropouts
