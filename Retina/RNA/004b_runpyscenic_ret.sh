#!/bin/bash

#PBS -N scenic
#PBS -l nodes=1:ppn=7
#PBS -l walltime=15:10:0
#PBS -l mem=20GB


#loading conda environment for pyscenic
source activate pyscn

##using loom file generated from prepforSCENIC to get adjcancey matrcex
pyscenic grn ${SCE}_ret.loom hg_tfs.txt -m genie3 -o adj${SCE}.csv --num_workers 8 --seed 10
#add correlation
pyscenic add_cor adj${SCE}.csv ${SCE}_ret.loom \
    --output cor${SCE}.tsv \
    --mask_dropouts

