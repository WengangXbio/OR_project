
GTreal1=$1
GTreal2=$2
GTobsv1=$3
GTobsv2=$4
vcfid=$5
GTrealy1=$(awk -F'\t' -v a=$GTreal1 -v b=$vcfid '$1==b && $2==a' all_bubble.ORgt |cut -f3-)
GTrealy2=$(awk -F'\t' -v a=$GTreal2 -v b=$vcfid '$1==b && $2==a' all_bubble.ORgt |cut -f3-)
GTobsvy1=$(awk -F'\t' -v a=$GTobsv1 -v b=$vcfid '$1==b && $2==a' all_bubble.ORgt |cut -f3-)
GTobsvy2=$(awk -F'\t' -v a=$GTobsv2 -v b=$vcfid '$1==b && $2==a' all_bubble.ORgt |cut -f3-)
GTrealy=$(echo -e "$GTrealy1,$GTrealy2")
GTobsvy=$(echo -e "$GTobsvy1,$GTobsvy2")
GTrealy=$(echo $GTrealy|tr ',' '\n' |grep -v "*" |tr '\n' ',' |sed 's/,$//')
GTobsvy=$(echo $GTobsvy|tr ',' '\n' |grep -v "*" |tr '\n' ',' |sed 's/,$//')

if [[ -n "$GTrealy" && -z "$GTobsvy" ]]; then
    CR=0
    FD="NA"
elif [[ -z "$GTrealy" && -n "$GTobsvy" ]]; then
    CR="NA"
    FD=1
elif [[ -z "$GTrealy" && -z "$GTobsvy" ]]; then
    CR="1"
    FD="0"
else
GTrealf=$(echo $GTrealy |tr ',' '\n' |sort |uniq -c |while IFS= read -r line; do
count=$(echo $line |awk '{print $1}')
type=$(echo $line |awk '{print $2}')
for n in `seq 1 $count` ; do
echo $type.$n
done
done)
GTobsvf=$(echo $GTobsvy |tr ',' '\n' |sort |uniq -c |while IFS= read -r line; do
count=$(echo $line |awk '{print $1}')
type=$(echo $line |awk '{print $2}')
for n in `seq 1 $count` ; do
echo $type.$n
done
done)
overlap=$(echo $GTrealf $GTobsvf |tr ' ' '\n' |sort |uniq -d |wc -l)
truenum=$(echo $GTrealf |tr ' ' '\n' |wc -l)
obsvnum=$(echo $GTobsvf |tr ' ' '\n' |wc -l)
CR=$(echo "scale=4; $overlap / $truenum" | bc)
FD=$(echo "scale=4; ($obsvnum - $overlap) / $obsvnum" | bc)
fi
echo -e "$vcfid\t$CR\t$FD"
