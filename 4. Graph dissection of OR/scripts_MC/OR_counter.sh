ORlist=$1
genomelist=$2
collapse=$3
codingbed=$4
batch=$5
grep -f ${ORlist} ${collapse} > test.collapse
awk 'BEGIN{OFS="\t"} {print $1,$3,$4,$2}' test.collapse | \
bedtools intersect -a ${codingbed} -b - -wa -wb > intersect
cat ${ORlist} > or.counts_${batch}.matrix
for genome0 in `cat ${genomelist}` ; do
grep "^${genome0}" intersect > test.intersect
score=0
for OR in `cat ${ORlist} ` ; do
matchOR=$(awk -v a=${OR} '$8==a {print $4}' test.intersect)
for gor in `echo $matchOR` ; do
awk -v a=$gor '$4==a' test.intersect |wc -l
done > xs1
while read -r number; do
    score=$(awk "BEGIN {print $score + 1/$number}")
done < xs1
echo -e $score
score=0
done |paste or.counts_${batch}.matrix - > temp
mv temp or.counts_${batch}.matrix
done 
cut -f4 intersect |sort |uniq > target.gene
rm xs1 test.intersect intersect test.collapse

#cut -f4 intersect |sort |uniq > target.gene
