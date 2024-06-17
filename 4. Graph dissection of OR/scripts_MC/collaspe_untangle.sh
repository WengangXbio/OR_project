inputprj=$1
inputseq=$2
outputfile=$3
for query in `cut -f1 ${inputseq}` ; do
awk -v a=$query '$1==a' ${inputprj} > gggenes
for subj in `cut -f2 gggenes |sort |uniq ` ; do
awk -v a=${subj} '$2==a' gggenes |cut -f3,4 |sed "s/\t/\n/g" |sort -k1,1n > strand_sort
start=$(head -n1 strand_sort)
t1=$(awk -v a=$start '$1 <= a+1200'  strand_sort |head -n1)
t2=$(awk -v a=$start '$1 <= a+1200'  strand_sort |tail -n1)
remain=$(awk -v a=$start '$1 >= a+1200' strand_sort)
while [ -n "$remain" ]; do 
echo $remain |sed "s/ /\n/g" > strand_sort
start=$(head -n1 strand_sort)
addt1=$(awk -v a=$start '$1 <= a+1200'  strand_sort |head -n1)
addt2=$(awk -v a=$start '$1 <= a+1200'  strand_sort |tail -n1)
t1="$t1 $addt1"
t2="$t2 $addt2"
remain=$(awk -v a=$start '$1 >= a+1200' strand_sort)
done
list1=($t1)
list2=($t2)
for i in "${!list1[@]}"; do
    echo -e "$query\t$subj\t${list1[i]}\t${list2[i]}"
done
done
done | awk '$4-$3 > 500' > ${outputfile}
rm strand_sort gggenes
#cut -f4 intersect |sort |uniq > target.gene

