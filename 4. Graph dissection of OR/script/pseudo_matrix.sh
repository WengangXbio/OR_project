ORlist=$1
genomelist=$2
collapse=$3
codingbed=$4
psedubed=$5
batch=$6
grep -f ${ORlist} ${collapse} > test.collapse
awk 'BEGIN{OFS="\t"} {print $1,$3,$4,$2}' test.collapse | \
bedtools intersect -a ${codingbed} -b - -wa -wb > intersect
awk 'BEGIN{OFS="\t"} {print $1,$3,$4,$2}' test.collapse | \
bedtools intersect -a ${psedubed} -b - -wa -wb > intersect.dead
cat ${ORlist} > or.pseudo_${batch}.matrix
for genome0 in `cat ${genomelist}` ; do
grep "^${genome0}" intersect > test.intersect
grep "^${genome0}" intersect.dead > test.intersect.dead
for OR in `cat ${ORlist} ` ; do
matchOR_coding=$(awk -v a=${OR} '$8==a' test.intersect |wc -l)
matchOR_dead=$(awk -v a=${OR} '$8==a' test.intersect.dead|wc -l)
if [ "$matchOR_coding" -gt 0 ]
then
echo 1
else
if [ "$matchOR_dead" -gt 0 ]
then
echo 0
else
echo -1
fi
fi
done |paste or.pseudo_${batch}.matrix - > temp
mv temp or.pseudo_${batch}.matrix
done

