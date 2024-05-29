fna=$1
echo "doing ${fna}"
awk '$7 >= 50000 && $10 >= 30000 && $12 == 60' ${fna}.paf |grep "^NC_" >  ${fna}.paf.gt10k
echo "subject filtering"
for cotig in ` awk '{print $6}' ${fna}.paf.gt10k |sort |uniq` ; do
mapn=$(awk -v a=${cotig} '$6==a' ${fna}.paf.gt10k |awk '{print $1}' |uniq |wc -l)
mapq=$(awk -v a=${cotig} '$6==a' ${fna}.paf.gt10k |awk '{print $1}' |uniq)
for mapqq in ` echo $mapq ` ; do
len=$(awk -v a=${cotig} '$6==a' ${fna}.paf.gt10k | awk -v a=${mapqq} '$1==a' | awk '{sum += $10} END {print sum}')
echo -e "$mapqq\t$len"
done > temp.len
maxq=$(sort -k2,2n temp.len |tail -n1|cut -f1)
awk -v a=${cotig} '$6==a' ${fna}.paf.gt10k | awk -v a=${maxq} '$1==a' |awk 'BEGIN{OFS="\t"} {print $6,$8,$9}' | sort -k2,2n > temp.max.bed
awk -v a=${cotig} '$6==a' ${fna}.paf.gt10k | awk -v a=${maxq} '$1!=a' |while IFS= read -r line; do
overlap=$(echo $line |awk 'BEGIN{OFS="\t"} {print $6,$8,$9,$1,$3,$4}' |bedtools intersect -a temp.max.bed -b - -wb |wc -l)
echo -e "$overlap\t$line"
done > temp.other.overlap
if [ "$mapn" = 1 ]
then
awk -v a=${cotig} '$6==a' ${fna}.paf.gt10k 
else
awk -F'\t' '$1==0' temp.other.overlap | cut -f2-
awk -v a=${cotig} '$6==a' ${fna}.paf.gt10k | awk -v a=${maxq} '$1==a'
fi
done > ${fna}.paf.gt10k.subject_filter

echo "query filtering"
for subj in ` awk '{print $1}' ${fna}.paf.gt10k.subject_filter |sort |uniq` ; do
mapn=$(grep $subj ${fna}.paf.gt10k.subject_filter |awk '{print $6}' |uniq |wc -l)
mapq=$(grep $subj ${fna}.paf.gt10k.subject_filter |awk '{print $6}' |uniq)
for mapqq in ` echo $mapq ` ; do
len=$(grep $subj ${fna}.paf.gt10k.subject_filter | awk -v a=${mapqq} '$6==a' | awk '{sum += $10} END {print sum}')
echo -e "$mapqq\t$len"
done > temp.len
maxq=$(sort -k2,2n temp.len |tail -n1|cut -f1)
grep $subj ${fna}.paf.gt10k.subject_filter | awk -v a=${maxq} '$6==a'  |awk 'BEGIN{OFS="\t"} {print $1,$3,$4}' | sort -k2,2n > temp.max.bed
grep $subj ${fna}.paf.gt10k.subject_filter | awk -v a=${maxq} '$6!=a'  |while IFS= read -r line; do
overlap=$(echo $line |awk 'BEGIN{OFS="\t"} {print $1,$3,$4,$6,$8,$9}' |bedtools intersect -a temp.max.bed -b - -wb |wc -l)
echo -e "$overlap\t$line"
done > temp.other.overlap
if [ "$mapn" = 1 ]
then
grep $subj ${fna}.paf.gt10k.subject_filter
else
awk -F'\t' '$1==0' temp.other.overlap | cut -f2-
grep $subj ${fna}.paf.gt10k.subject_filter |awk -v a=${maxq} '$6==a' 
fi
done > ${fna}.paf.gt10k.subject_filter.query_filter

echo "build alignment map"
for cotig in ` awk '{print $6}' ${fna}.paf.gt10k.subject_filter.query_filter |sort |uniq` ; do
mapn=$(awk -v a=${cotig} '$6==a' ${fna}.paf.gt10k.subject_filter.query_filter |awk '{print $1}' |uniq |wc -l)
mapq=$(awk -v a=${cotig} '$6==a' ${fna}.paf.gt10k.subject_filter.query_filter |awk '{print $1}' |uniq)
multi=multi
if [ "$mapn" = 1 ]
then
echo -e "$fna\t$cotig\t$mapq"
else
echo -e "$fna\t$cotig\t$multi"
fi 
done > align.list

for cotig in `grep "multi" align.list |cut -f2 |sort |uniq ` ; do 
awk -v a=${cotig} '$6==a' ${fna}.paf.gt10k.subject_filter.query_filter |awk -v a=$fna 'BEGIN{OFS="\t"} {print $6,$8,$9,$1} '
done > align_multi.bed

for seq in ` awk '{print $4}' align_multi.bed |sort |uniq` ; do
grep $seq align_multi.bed |sort -k1,1 -k2,2n |bedtools merge -i - |awk -v a=$fna -v b=$seq 'BEGIN{OFS="\t"} {print a,$1,$2,$3,b} '
done > kmap1
Rscript kmap.merge.r
grep -v "multi" align.list |awk 'BEGIN{OFS="\t"} {print $1,$2,0,0,$3}' > kmap2
cat kmap1.merge kmap2 > $fna.align_map

