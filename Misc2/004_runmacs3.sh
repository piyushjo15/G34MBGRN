#!/bin/bash
## This script calls peak summits for TF-bindings per sample-TF combination using macs3
module load SAMtools

source /path/to/macs3/bin/activate

TMP=/path/to/TMP
OUT=/path/to/output

##bwa
BAM=$OUT/${sample}_${tf}
INPUT=$OUT/MBinput_Input/MBinput_Input.sorted.dedup.bam ##pooled ig input control

macs3 callpeak -t $BAM/${sample}_${tf}.sorted.dedup.bam -c $INPUT \
-f BAMPE -g hs -p 0.001 --tempdir $TMP --seed 10 --outdir $OUT/macs3/ -n ${sample}_${tf}
