#!/bin/bash

## This script is for alignment of ChIP-seq FASTQ files.
## BWA aligner is used as sequence length is 50 bp

module load BWA/
module load SAMtools
module load sambamba


FQS=/path/to/fastqs/
REF=/path/to/bwa_ref/
TMP=/path/to/TMP
OUT=/path/to/output

cd $OUT

mkdir -p $OUT/${sample}_${tf}
echo "bwa alignment for ..."
echo ${sample}_${tf} ## sample and TF combination

cd $IN/${sample}_${tf}
##alignment
echo "performing alignment..."
bwa aln -t 6 $REF  $FQS/${sample}_${tf}/*_R1.fastq.gz >  ${sample}_${tf}_1.sai
bwa aln -t 6 $REF  $FQS/${sample}_${tf}/*_R2.fastq.gz >  ${sample}_${tf}_2.sai

##same
bwa sampe $REF ${sample}_${tf}_1.sai  ${sample}_${tf}_2.sai $FQS/*_R1.fastq.gz $FQS/*_R2.fastq.gz > ${sample}_${tf}.sam

##sam to bam
samtools view -h -S -b -o ${sample}_${tf}.bam ${sample}_${tf}.sam

rm ${sample}_${tf}.sam

## sorting
sambamba sort -t 6 -o ${sample}_${tf}.sorted.bam ${sample}_${tf}.bam --tmpdir=$TMP

## mark duplicates
sambamba markdup -t 6 ${sample}_${tf}.sorted.bam ${sample}_${tf}.sorted.mdup.bam --tmpdir=$TMP

## remove multimapped, unmapped and duplicated
sambamba view -h -t 6 -f bam \
	-F "[XS] == null and not unmapped  and not duplicate" \
	${sample}_${tf}.sorted.mdup.bam > ${sample}_${tf}.sorted.dedup.bam


## indexing deduped bam
samtools index ${sample}_${tf}.sorted.dedup.bam



