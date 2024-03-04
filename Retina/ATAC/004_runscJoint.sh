#!/bin/bash

#BSUB -J Run_scJoint
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 1:00
#BSUB -n 1
#BSUB -R "rusage[mem=30GB]"

cd $HOME/MBsnANA/Retina/ATAC/ArchR/scJoint
module load miniconda/4.9.2
source activate pjexv4
echo "Running scJoint ..."
python main_scjoint.py ##this requires config.py which is Retina-config.py
python process_scJoint.py

