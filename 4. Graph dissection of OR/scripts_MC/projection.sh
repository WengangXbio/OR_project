#BSUB -J OMGproj
#BSUB -n 10
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=4GB]
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -q normal


source activate ~/.conda/envs/odgi
NC=OMG
og=cattle31.${NC}.full.og
odgi inject -i ${og} -b ${NC}.func.proc.bed -o injref.og -t 10
cut -f4 ${NC}.func.proc.bed > coding.name
odgi untangle -R coding.name -i injref.og -g -j 0.1 -n 100 -t 10 > injref.untangle.tsv

