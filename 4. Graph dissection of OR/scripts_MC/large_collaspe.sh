#BSUB -J blast
#BSUB -n 1
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=4GB]
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -q normal


module load anaconda3/4.12.0
source activate odgi 
./collaspe_untangle.sh injref.untangle.tsv inputseq injref.untangle.tsv.collapse

