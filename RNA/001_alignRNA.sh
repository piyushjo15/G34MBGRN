#!/bin/bash


module load samtools
export PATH=$PATH:$HOME/STAR/source


##define input and output directories
IDX=$HOME/path/to/starindex #index for star aligner
OUT=$HOME/path/to/STARoutput #output directory
WL=$HOME/CR_ARC/lib/python/cellranger/barcodes/737K-arc-v1.txt ##cellranger barcodes
FQ=$HOME/path/to/${SCE} ##path to FASTQs, SCE= sample id
cd $FQ
##C1 are all R1 fastqs and C2 are all R2 fastqs
C1=`ls -m *${SCE}*_R1_001.fastq.gz | tr -d '\n' | tr -d ' '`
C2=`ls -m *${SCE}*_R2_001.fastq.gz | tr -d '\n' | tr -d ' '`

##Run STARsolo
STAR --runMode alignReads --soloType CB_UMI_Simple \
--quantMode GeneCounts --runThreadN 8 \
--soloUMIlen 12 --soloCBwhitelist $WL \
--readFilesIn $C2 $C1 \
--genomeDir $IDX \
--twopassMode Basic \
--outFileNamePrefix $OUT/${SCE}/${SCE} \
--readFilesCommand zcat \
--soloFeatures Gene GeneFull \
--soloUMIfiltering MultiGeneUMI \
--soloCBmatchWLtype 1MM_multi_pseudocounts \
--outSAMtype BAM Unsorted \
--soloCellFilter None \
--outSAMunmapped Within \
--outSAMmultNmax 1 \
--limitSjdbInsertNsj 1500000
