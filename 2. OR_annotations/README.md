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

### Evaluation of OR regions assembly by distribution and counts
```
for genome in `ls *.fna`  ; do
numb=$(cut -f1 ${genome}.g2or_itera2_func.bed | uniq |wc -l)
numa=$(cat ${genome}.g2or_itera2_func.bed |wc -l)
echo -e "$numa\t$numb\t$genome"
done
```

| Genome                                                             | OR counts | Distribution (by chromosome/ scaffold) |
|--------------------------------------------------------------------|-----------|----------------------------------------|
| GCA_002933975.1_ASM293397v1_genomic.fna                            | 1028      | 23                                     |
| GCA_003369685.2_UOA_Angus_1_genomic.fna                            | 1098      | 28                                     |
| GCA_009493655.1_ARS_UNL_Btau-highland_paternal_1.0_alt_genomic.fna | 1097      | 22                                     |
| GCA_021234555.1_ARS-LIC_NZ_Jersey_genomic.fna                      | 1071      | 26                                     |
| GCA_021347905.1_ARS-LIC_NZ_Holstein-Friesian_1_genomic.fna         | 1034      | 41                                     |
| GCA_028973685.2_SNU_Hanwoo_2.0_genomic.fna                         | 1113      | 24                                     |
| GCA_029378745.1_NIAB-ARS_B.indTharparkar_mat_pri_1.0_genomic.fna   | 1129      | 21                                     |
| GCA_030267345.1_ASM3026734v1_genomic.fna                           | 1106      | 21                                     |
| GCA_030267355.1_ASM3026735v1_genomic.fna                           | 1121      | 23                                     |
| GCA_030267375.1_ASM3026737v1_genomic.fna                           | 1104      | 23                                     |
| GCA_030269505.1_ASM3026950v1_genomic.fna                           | 1105      | 22                                     |
| GCA_030269815.1_ASM3026981v1_genomic.fna                           | 1099      | 22                                     |
| GCA_030270685.1_ASM3027068v1_genomic.fna                           | 1109      | 21                                     |
| GCA_030270715.1_ASM3027071v1_genomic.fna                           | 1116      | 21                                     |
| GCA_030271795.1_ASM3027179v1_genomic.fna                           | 1132      | 22                                     |
| GCA_030271805.1_ASM3027180v1_genomic.fna                           | 1094      | 22                                     |
| GCA_030272135.1_ASM3027213v1_genomic.fna                           | 1110      | 23                                     |
| GCA_034097375.1_YAU_Btau_1.0_genomic.fna                           | 1142      | 23                                     |
| GCA_905123515.1_ROSLIN_BTT_NDA1_genomic.fna                        | 1043      | 25                                     |
| GCA_905123885.1_ROSLIN_BTI_ANK1_genomic.fna                        | 1076      | 117                                    |
| GCA_947034695.1_seqoccin.Bt.char.v1.0_genomic.fna                  | 1100      | 22                                     |
| GCF_000247795.1_Bos_indicus_1.0_genomic.fna                        | 171       | 15                                     |
| GCF_002263795.3_ARS-UCD2.0_genomic.fna                             | 1079      | 25                                     |
| GCF_003369695.1_UOA_Brahman_1_genomic.fna                          | 1128      | 27                                     |
| GCA_005887515.3_BosGru3.1_genomic.fna                              | 697       | 24                                     |
| GCA_007844835.1_NRC_Mithun_1_genomic.fna                           | 768       | 240                                    |
| GCA_009493645.1_ARS_UNL_BGru_maternal_1.0_p_genomic.fna            | 1058      | 25                                     |
| GCA_014182915.2_ARS_UOA_Gaur_1.1_genomic.fna                       | 1048      | 43                                     |
| GCA_027580195.1_NWIPB_WYAK_1.0_genomic.fna                         | 832       | 36                                     |
| GCA_027580245.1_NWIPB_DYAK_1.0_genomic.fna                         | 844       | 38                                     |
| GCA_946052875.1_ARS_Pied_mat1.0_genomic.fna                        | 1127      | 21                                     |
| GCA_963879515.1_ETH_BisBon1_genomic.fna                            | 1100      | 23                                     |
| GCF_000754665.1_Bison_UMD1.0_genomic.fna                           | 877       | 311                                    |
| GCF_032452875.1_ARS-OSU_banteng_1.0_genomic.fna                    | 1147      | 25                                     |
