#!/bin/bash

##script for running homer motif enrichment

cd /path/to/output/macs3/

OUT=motifenr
TMP=/path/to/TMP

module load homer/4.11

findMotifsGenome.pl EOMES.com.bed hg38r $OUT/EOMES -p 8 -preparsedDir $TMP -fdr 
findMotifsGenome.pl CRX.com.bed hg38r $OUT/CRX -p 8 -preparsedDir $TMP -fdr 
findMotifsGenome.pl PAX6.bed hg38r $OUT/PAX6 -p 8 -preparsedDir $TMP -fdr 
