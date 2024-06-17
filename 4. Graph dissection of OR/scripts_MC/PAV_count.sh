ORlist=$1
genomelist=$2
collapse=$3
codingbed=$4
batch=$5
grep -f ${ORlist} ${collapse} > test.collapse
awk 'BEGIN{OFS="\t"} {print $1,$3,$4,$2}' test.collapse | \
bedtools intersect -a ${codingbed} -b - -wa -wb > intersect

cat ${ORlist} > or.PAV_${batch}.matrix
for genome0 in `cat ${genomelist}` ; do
for OR in `cat ${ORlist} ` ; do
grep "^${genome0}" injref.untangle.tsv.collapse |awk -v a=${OR} '$2==a' |wc -l
done |paste or.PAV_${batch}.matrix - > temp
mv temp or.PAV_${batch}.matrix
done
cut -f4 intersect |sort |uniq > target.gene

