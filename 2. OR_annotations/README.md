## OR annotation by genome2OR
```
cd /share/home/yzwl_zhangwg/OR_project/genome2OR
for genome in `ls *.fna`  ; do
outdir="${genome}.g2or"
sed "s/INPUT/${genome}/g" genome2OR.sh | sed "s/OUTPUT/${outdir}/g" > ${genome}.sh
csub < ${genome}.sh
done
```
### Convert genome2OR results into bed
```
for genome in `ls *.fna`  ; do
samtools faidx ${genome}.g2or/ORannotation_itera2_final_func_dna_ORs.fasta 
awk -F'_' 'BEGIN{OFS="\t"} {print $2,$3,$4,$5}' ${genome}.g2or/ORannotation_itera2_final_func_dna_ORs.fasta.fai |sort -k1,1 -k2,2n > ${genome}.g2or/ORannotation_itera2_final_func_dna_ORs.bed
sh modify_coor.sh ${genome}.g2or/ORannotation_itera2_final_func_dna_ORs.bed ${genome}.g2or/ORannotation_itera2_final_func_dna_ORs.modify.bed
sed "s/@#@/_/g" ${genome}.g2or/ORannotation_itera2_final_func_dna_ORs.modify.bed > ${genome}.g2or/ORannotation_GCF_002263795.3_interation.bed
cp ${genome}.g2or/ORannotation_GCF_002263795.3_interation.bed ${genome}.g2or_itera2_func.bed
done
```
```
[yzwl_zhangwg@mgt15 genome2OR]$ cat genome2OR.sh 
#!/bin/bash
#CSUB -J genome2OR
#CSUB -q cpu
#CSUB -o %J.out
#CSUB -e %J.error
#CSUB -n 24   #core
#CSUB -R span[hosts=1]

module load anaconda3/4.12.0
source activate genome2or
Iteration Mammalia OUTPUT INPUT 
```
