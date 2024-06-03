#!/bin/bash
#CSUB -J odgiinj
#CSUB -q cpu
#CSUB -o %J.out
#CSUB -e %J.error
#CSUB -n 12   #core
#CSUB -R span[hosts=1]

module load anaconda3/4.12.0
source activate odgi 
og=OMG
odgi inject -i ${og} -b coding_annotation.bed -o injref.og -t 24
cut -f4 coding_annotation.bed > coding.name
odgi untangle -R coding.name -i injref.og -g -j 0.1 -n 100 -t 24 > injref.untangle.tsv

