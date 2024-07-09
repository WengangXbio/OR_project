cut -f1 combined_ORgt.txt |sort |uniq | while IFS= read -r orid; do
awk -v a=$orid '$1==a {print $3}' combined_ORgt.txt |grep "OR"|tr ',' '\n' |sort |uniq | awk -v a=$orid '{print a"_OR"NR,$0}' > temp
for i in $(seq 1 100); do
pair=$(awk '{print $1}' temp |sort |uniq |shuf |head -2 |sort )
if [ $(awk '{print $1}' temp |sort |uniq |wc -l) -gt 1 ]; then
gt1=$(echo $pair |awk '{print $1}')
gt2=$(echo $pair |awk '{print $2}')
gtor1=$(awk -v a=$gt1 '$1==a {print $2}' temp  |cut -c 2- |tr '|' '\n' |grep OR|sort |uniq)
gtor2=$(awk -v a=$gt2 '$1==a {print $2}' temp  |cut -c 2- |tr '|' '\n' |grep OR|sort |uniq)
same=$(echo $gtor1 $gtor2 |tr ' ' '\n' |sort |uniq -c |awk '$1>1' |wc -l)
if [ "$same" -gt 0 ]; then
awk -v a="$gt2" -v b="$gt1" 'BEGIN{FS=OFS=" "}{if ($1==a) $1=b; print}' temp > temp0
mv temp0 temp
fi
fi
done
cat temp
done > combined_ORgt.ORID.track

cut -f1 combined_ORgt.txt |sort |uniq | while IFS= read -r orid; do
awk -v a=$orid '$1==a' combined_ORgt.txt |grep "OR" |while IFS= read -r line; do
gt=$(echo $line |awk '{print $2}')
orlt=$(echo $line |awk '{print $3}')
for or in `echo $orlt |tr ',' '\n'` ; do
grep $orid combined_ORgt.ORID.track |awk -v a=$or '$2==a {print $1}'
done |sort |uniq -c |awk -v a=$gt -v b=$orid '{print a,b,$2,$1}'
done
done > combined_ORgt.dosage
