#!/bin/bash
#CSUB -J odgiinj
#CSUB -q cpu
#CSUB -o %J.out
#CSUB -e %J.error
#CSUB -n 2   #core
#CSUB -R span[hosts=1]

module load anaconda3/4.12.0
source activate odgi 
./collaspe_untangle.sh injref.untangle.tsv inputseq injref.untangle.tsv.collapse

