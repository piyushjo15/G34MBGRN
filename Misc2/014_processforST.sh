
#!/bin/bash

## This sciprt uses hic-pro output to produce output files for SpectralTAD input
## and obtain interaction loops
## It requires HiCLift, cooler and peakachu
cd HIC/GSE240410_data
  
echo "processing sample for SpectralTAD input .."
chromsizes= /path/to/chromsize.bed
## convert to .mcool format 
HiCLift --input MB288.allValidPairs.txt.gz --input-format hic-pro --output-format cool \
--out-pre MB288_HL --out-chromsizes ${chromsizes} --in-assembly hg38 --out-assembly hg38 --memory 40G

## obtain input for SpectralTAD
cooler dump --join MB288_HL.mcool::resolutions/10000 > MB288_HL.10kb.txt
cooler dump --join MB288_HL.mcool::resolutions/25000 > MB288_HL.25kb.txt

cd ..
mkdir CLOOPs
## run peackahu
##25kb
peakachu score_genome -r 25000 --clr-weight-name weight \
-p GSE240410_data/MB288_HL.mcool::resolutions/25000 \
-O CLOOPS/MB288_25kb_scores.bedpe -m high-confidence.400million.25kb.w5.pkl
peakachu pool -r 25000 -i CLOOPS/MB288_25kb_scores.bedpe -o CLOOPS/MB288_25kb.loops.0.95.bedpe -t 0.949

##10kb
peakachu score_genome -r 10000 --clr-weight-name weight \
-p GSE240410_data/MB288_HL.mcool::resolutions/10000 \
-O CLOOPS/MB288_10kb_scores.bedpe -m high-confidence.400million.10kb.w6.pkl
peakachu pool -r 10000 -i CLOOPS/MB288_10kb_scores.bedpe -o CLOOPS/MB288_10kb.loops.0.95.bedpe -t 0.949
