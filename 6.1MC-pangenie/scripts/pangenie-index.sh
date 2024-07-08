#!/bin/bash
#CSUB -J pangenie-index
#CSUB -q cpu
#CSUB -o %J.out
#CSUB -e %J.error
#CSUB -n 24   #core
#CSUB -R span[hosts=1]

module load anaconda3/4.12.0
source activate pangenie 
/share/home/yzwl_zhangwg/tools/pangenie/build/src/PanGenie-index -r GCF_002263795.3_ARS-UCD2.0_genomic.fna -v cattle31.raw.vcfbub.vcf.diploid.noheader.removeN.vcf -o preprocessing -e 100000 -t 16
