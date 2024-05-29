#!/bin/bash
#CSUB -J CHRNkmap
#CSUB -q cpu
#CSUB -o %J.out
#CSUB -e %J.error
#CSUB -n 6   #core
#CSUB -R span[hosts=1]


module load anaconda3/4.12.0
source activate minimap2

wd="/share/home/yzwl_zhangwg/OR_project/PGGB/pan_by_chr_new_new"
chrn=CHRN
rm $wd/$chrn -rf
mkdir $wd/$chrn -p
ass=$(awk -v a=$chrn '$1==a {print $2}' ARS20_chr_NC.info)

./kmap.step2.sh $chrn $ass
