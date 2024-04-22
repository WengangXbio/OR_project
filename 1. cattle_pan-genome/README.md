# 1. Cattle Pan-genome
## 1.1 Assembly summary

| shortnames          | Genome                                                             | Linages                        |
|---------------------|--------------------------------------------------------------------|--------------------------------|
| UOA_Angus_1         | GCA_003369685.2_UOA_Angus_1_genomic.fna                            | Bos taurus                     |
| Btau_highland       | GCA_009493655.1_ARS_UNL_Btau-highland_paternal_1.0_alt_genomic.fna | Bos taurus                     |
| Jersey              | GCA_021234555.1_ARS-LIC_NZ_Jersey_genomic.fna                      | Bos taurus                     |
| Holstein_Friesian   | GCA_021347905.1_ARS-LIC_NZ_Holstein-Friesian_1_genomic.fna         | Bos taurus                     |
| Hanwoo              | GCA_028973685.2_SNU_Hanwoo_2.0_genomic.fna                         | Bos taurus                     |
| ARS_B_indTharparkar | GCA_029378745.1_NIAB-ARS_B.indTharparkar_mat_pri_1.0_genomic.fna   | Bos indicus                    |
| Weining             | GCA_030267345.1_ASM3026734v1_genomic.fna                           | Bos indicus                    |
| Xiangxi             | GCA_030267355.1_ASM3026735v1_genomic.fna                           | Bos indicus                    |
| Yiling              | GCA_030267375.1_ASM3026737v1_genomic.fna                           | Bos indicus                    |
| Jinjiang            | GCA_030269505.1_ASM3026950v1_genomic.fna                           | Bos indicus                    |
| Leiqiong            | GCA_030269815.1_ASM3026981v1_genomic.fna                           | Bos indicus                    |
| Dabieshan           | GCA_030270685.1_ASM3027068v1_genomic.fna                           | Bos indicus                    |
| Guanling            | GCA_030270715.1_ASM3027071v1_genomic.fna                           | Bos indicus                    |
| Wenshangaofeng      | GCA_030271795.1_ASM3027179v1_genomic.fna                           | Bos indicus                    |
| Weizhou             | GCA_030271805.1_ASM3027180v1_genomic.fna                           | Bos indicus                    |
| Lincanggaofeng      | GCA_030272135.1_ASM3027213v1_genomic.fna                           | Bos indicus                    |
| Nelore              | GCF_000247795.1_Bos_indicus_1.0_genomic.fna                        | Bos indicus                    |
| ROSLIN_BTI_ANK1     | GCA_905123885.1_ROSLIN_BTI_ANK1_genomic.fna                        | Bos taurus                     |
| ROSLIN_BTT_NDA1     | GCA_905123515.1_ROSLIN_BTT_NDA1_genomic.fna                        | Bos taurus                     |
| Yunling             | GCA_034097375.1_YAU_Btau_1.0_genomic.fna                           | Bos indicus                    |
| UOA_Brahman_1       | GCF_003369695.1_UOA_Brahman_1_genomic.fna                          | Bos taurus                     |
| ARS_UCD2_0          | GCF_002263795.3_ARS-UCD2.0_genomic.fna                             | Bos taurus                     |
| Gir                 | GCA_002933975.1_ASM293397v1_genomic.fna                            | Bos indicus                    |
| BosGru3_1           | GCA_005887515.3_BosGru3.1_genomic.fna                              | Bos grunniens (yak)            |
| Gayal               | GCA_007844835.1_NRC_Mithun_1_genomic.fna                           | Bos frontalis (gayal)          |
| BGru_maternal       | GCA_009493645.1_ARS_UNL_BGru_maternal_1.0_p_genomic.fna            | Bos grunniens (yak)            |
| Gaur                | GCA_014182915.2_ARS_UOA_Gaur_1.1_genomic.fna                       | Bos gaurus (gaur)              |
| WYAK                | GCA_027580195.1_NWIPB_WYAK_1.0_genomic.fna                         | Bos mutus (wild yak)           |
| DYAK                | GCA_027580245.1_NWIPB_DYAK_1.0_genomic.fna                         | Bos grunniens (yak)            |
| ARS_Pied_mat1       | GCA_946052875.1_ARS_Pied_mat1.0_genomic.fna                        | Bos gaurus (gaur)              |
| seqoccin_Bt_char    | GCA_947034695.1_seqoccin.Bt.char.v1.0_genomic.fna                  | Bos taurus                     |
| European_bison      | GCA_963879515.1_ETH_BisBon1_genomic.fna                            | Bison bonasus (European bison) |
| American_bison      | GCF_000754665.1_Bison_UMD1.0_genomic.fna                           | Bison bison bison              |
| Banteng             | GCF_032452875.1_ARS-OSU_banteng_1.0_genomic.fna                    | Bos javanicus (banteng)        |

## 1.2 Draw distribution of OR genes in cattle and human
### 1.2.1 Extract OR gene gff (cattle and human)
```
#cattle
grep "olfactory receptor" GCF_002263795.3_ARS-UCD2.0_genomic.gff |awk -F'\t' '$3=="gene"||$3=="pseudogene"' > GCF_002263795.3_ARS-UCD2.0_genomic.OR.body
grep "^#" GCF_002263795.3_ARS-UCD2.0_genomic.gff > GCF_002263795.3_ARS-UCD2.0_genomic.OR.header
cat GCF_002263795.3_ARS-UCD2.0_genomic.OR.header GCF_002263795.3_ARS-UCD2.0_genomic.OR.body > GCF_002263795.3_ARS-UCD2.0_genomic.OR.gff
#human
grep "olfactory receptor" GCF_000001405.40_GRCh38.p14_genomic.gff |awk -F'\t' '$3=="gene"||$3=="pseudogene"' > GCF_000001405.40_GRCh38.p14_genomic.OR.body
grep "^#" GCF_000001405.40_GRCh38.p14_genomic.OR.body > GCF_000001405.40_GRCh38.p14_genomic.OR.header
cat GCF_000001405.40_GRCh38.p14_genomic.OR.header GCF_000001405.40_GRCh38.p14_genomic.OR.body > GCF_000001405.40_GRCh38.p14_genomic.OR.gff
```
### 1.2.2 [Plot in R](https://github.com/WengangXbio/OR_project/blob/f24859a6f5c74839a430af168e71fc48c5621a34/1.%20cattle_pan-genome/0.%20script/Plot_distribution.r) 
OR distribution in cattle ARS2.0
![cattle](https://github.com/WengangXbio/OR_project/blob/b47049e7e42c21cef1a6b1acb78a784550ff9006/1.%20cattle_pan-genome/1.%20Figure/OR_distribution_cattle.png)
OR distribution in human GhRC38
![human](https://github.com/WengangXbio/OR_project/blob/b47049e7e42c21cef1a6b1acb78a784550ff9006/1.%20cattle_pan-genome/1.%20Figure/OR_distribution.png)

