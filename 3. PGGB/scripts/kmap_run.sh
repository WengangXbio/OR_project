#!/bin/bash
#CSUB -J ASSEMBLYkmap
#CSUB -q cpu
#CSUB -o %J.out
#CSUB -e %J.error
#CSUB -n 6   #core
#CSUB -R span[hosts=1]


module load anaconda3/4.12.0
source activate minimap2

fna="ASSEMBLY"
mkdir ${fna} -p
cp ./kmap.step1.sh ${fna}/kmap.step1.sh
cp ${fna}.paf ${fna}/${fna}.paf
cp kmap.merge.r ${fna}/kmap.merge.r
cd ${fna}
chmod +x ./kmap.step1.sh
sh ./kmap.step1.sh ${fna}
