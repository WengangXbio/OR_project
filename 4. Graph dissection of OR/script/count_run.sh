#!/bin/bash
#CSUB -J countCHRN
#CSUB -q cpu
#CSUB -o %J.out
#CSUB -e %J.error
#CSUB -n 2   #core
#CSUB -R span[hosts=1]

module load anaconda3/4.12.0
source activate odgi 
nchr=CHRN
cp /share/home/yzwl_zhangwg/OR_project/PGGB/script/count.matrix.sh ./count.matrix.sh
chmod +x /count.matrix.sh
og=chr${nchr}.final.og
odgi paths -i ${og} -L > inputseq
./collaspe_untangle.sh injref.untangle.tsv inputseq injref.untangle.tsv.collapse # collaspe_untangle.sh input.tsv input.seq output.tsv
#./count.matrix.sh ${nchr}
