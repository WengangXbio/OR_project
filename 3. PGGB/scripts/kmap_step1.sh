#Build consolidate alignment map corresponding to reference for each assembly
fnalist="GCA_002933975.1_ASM293397v1_genomic.fna GCA_003369685.2_UOA_Angus_1_genomic.fna GCA_005887515.3_BosGru3.1_genomic.fna GCA_007844835.1_NRC_Mithun_1_genomic.fna GCA_009493645.1_ARS_UNL_BGru_maternal_1.0_p_genomic.fna GCA_009493655.1_ARS_UNL_Btau-highland_paternal_1.0_alt_genomic.fna GCA_014182915.2_ARS_UOA_Gaur_1.1_genomic.fna GCA_021234555.1_ARS-LIC_NZ_Jersey_genomic.fna GCA_021347905.1_ARS-LIC_NZ_Holstein-Friesian_1_genomic.fna GCA_027580195.1_NWIPB_WYAK_1.0_genomic.fna GCA_027580245.1_NWIPB_DYAK_1.0_genomic.fna GCA_028973685.2_SNU_Hanwoo_2.0_genomic.fna GCA_029378745.1_NIAB-ARS_B.indTharparkar_mat_pri_1.0_genomic.fna GCA_030267345.1_ASM3026734v1_genomic.fna GCA_030267355.1_ASM3026735v1_genomic.fna GCA_030267375.1_ASM3026737v1_genomic.fna GCA_030269505.1_ASM3026950v1_genomic.fna GCA_030269815.1_ASM3026981v1_genomic.fna GCA_030270685.1_ASM3027068v1_genomic.fna GCA_030270715.1_ASM3027071v1_genomic.fna GCA_030271795.1_ASM3027179v1_genomic.fna GCA_030271805.1_ASM3027180v1_genomic.fna GCA_030272135.1_ASM3027213v1_genomic.fna GCA_034097375.1_YAU_Btau_1.0_genomic.fna GCA_905123515.1_ROSLIN_BTT_NDA1_genomic.fna GCA_905123885.1_ROSLIN_BTI_ANK1_genomic.fna GCA_946052875.1_ARS_Pied_mat1.0_genomic.fna GCA_947034695.1_seqoccin.Bt.char.v1.0_genomic.fna GCA_963879515.1_ETH_BisBon1_genomic.fna GCF_000247795.1_Bos_indicus_1.0_genomic.fna GCF_000754665.1_Bison_UMD1.0_genomic.fna GCF_002263795.3_ARS-UCD2.0_genomic.fna GCF_003369695.1_UOA_Brahman_1_genomic.fna GCF_032452875.1_ARS-OSU_banteng_1.0_genomic.fna"
for fna in `echo $fnalist ` ; do
echo "doing ${fna}"
awk '$7 >= 50000 && $10 >= 30000 && $12 == 60' ${fna}.paf |grep "^NC_" >  ${fna}.paf.gt10k
for cotig in ` awk '{print $6}' ${fna}.paf.gt10k |sort |uniq` ; do
mapn=$(grep $cotig ${fna}.paf.gt10k |awk '{print $1}' |uniq |wc -l)
mapq=$(grep $cotig ${fna}.paf.gt10k |awk '{print $1}' |uniq)
for mapqq in ` echo $mapq ` ; do
len=$(grep $cotig ${fna}.paf.gt10k | grep $mapqq | awk '{sum += $10} END {print sum}')
echo -e "$mapqq\t$len"
done > temp.len
maxq=$(sort -k2,2n temp.len |tail -n1|cut -f1)
grep $cotig ${fna}.paf.gt10k | grep $maxq |awk 'BEGIN{OFS="\t"} {print $6,$8,$9}' | sort -k2,2n > temp.max.bed
grep $cotig ${fna}.paf.gt10k | grep -v $maxq |while IFS= read -r line; do
overlap=$(echo $line |awk 'BEGIN{OFS="\t"} {print $6,$8,$9,$1,$3,$4}' |bedtools intersect -a temp.max.bed -b - -wb |wc -l)
echo -e "$overlap\t$line"
done > temp.other.overlap
if [ "$mapn" = 1 ]
then
grep $cotig ${fna}.paf.gt10k
else
awk -F'\t' '$1==0' temp.other.overlap | cut -f2-
grep $cotig ${fna}.paf.gt10k |grep $maxq 
fi
done > ${fna}.paf.gt10k.subject_filter

for subj in ` awk '{print $1}' ${fna}.paf.gt10k.subject_filter |sort |uniq` ; do
mapn=$(grep $subj ${fna}.paf.gt10k.subject_filter |awk '{print $6}' |uniq |wc -l)
mapq=$(grep $subj ${fna}.paf.gt10k.subject_filter |awk '{print $6}' |uniq)
for mapqq in ` echo $mapq ` ; do
len=$(grep $subj ${fna}.paf.gt10k.subject_filter | grep $mapqq | awk '{sum += $10} END {print sum}')
echo -e "$mapqq\t$len"
done > temp.len
maxq=$(sort -k2,2n temp.len |tail -n1|cut -f1)
grep $subj ${fna}.paf.gt10k.subject_filter | grep $maxq |awk 'BEGIN{OFS="\t"} {print $1,$3,$4}' | sort -k2,2n > temp.max.bed
grep $subj ${fna}.paf.gt10k.subject_filter | grep -v $maxq |while IFS= read -r line; do
overlap=$(echo $line |awk 'BEGIN{OFS="\t"} {print $1,$3,$4,$6,$8,$9}' |bedtools intersect -a temp.max.bed -b - -wb |wc -l)
echo -e "$overlap\t$line"
done > temp.other.overlap
if [ "$mapn" = 1 ]
then
grep $subj ${fna}.paf.gt10k.subject_filter
else
awk -F'\t' '$1==0' temp.other.overlap | cut -f2-
grep $subj ${fna}.paf.gt10k.subject_filter |grep $maxq 
fi
done > ${fna}.paf.gt10k.subject_filter.query_filter

for cotig in ` awk '{print $6}' ${fna}.paf.gt10k.subject_filter.query_filter |sort |uniq` ; do
mapn=$(grep $cotig ${fna}.paf.gt10k.subject_filter.query_filter |awk '{print $1}' |uniq |wc -l)
mapq=$(grep $cotig ${fna}.paf.gt10k.subject_filter.query_filter |awk '{print $1}' |uniq)
multi=multi
if [ "$mapn" = 1 ]
then
echo -e "$fna\t$cotig\t$mapq"
else
echo -e "$fna\t$cotig\t$multi"
fi 
done > align.list

grep "multi" align.list |cut -f2 |sort |uniq | grep -f - ${fna}.paf.gt10k.subject_filter.query_filter |\
awk -v a=$fna 'BEGIN{OFS="\t"} {print $6,$8,$9,$1} ' > align_multi.bed
for seq in ` awk '{print $4}' align_multi.bed |sort |uniq` ; do
grep $seq align_multi.bed |sort -k1,1 -k2,2n |bedtools merge -i - |awk -v a=$fna -v b=$seq 'BEGIN{OFS="\t"} {print a,$1,$2,$3,b} '
done > kmap1
Rscript kmap.merge.r
grep -v "multi" align.list |awk 'BEGIN{OFS="\t"} {print $1,$2,0,0,$3}' > kmap2
cat kmap1.merge kmap2 > $fna.align_map
done
