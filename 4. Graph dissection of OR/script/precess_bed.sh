seqfile=$1
annofile=$2
cat $seqfile | while read line; do
#echo $line
genome1=$(echo ${line} | awk -F'[#]' '{print $1}')
anno=$annofile
seq=$(echo ${line} | awk -F'[#]' '{print $3}' |awk -F'[_]' '{print $1}' )
if [ "$seq" == NC ]; then
seq=$(echo ${line} | awk -F'[#]' '{print $3}' |awk -F'[_]' '{print $1"_"$2}' )
start=$(echo ${line} | awk -F'[#]' '{print $3}' |awk -F'[_]' '{print $3}')
end=$(echo ${line} | awk -F'[#]' '{print $3}' |awk -F'[_]' '{print $4}')
else
seq=$(echo ${line} | awk -F'[#]' '{print $3}' |awk -F'[_]' '{print $1}' )
start=$(echo ${line} | awk -F'[#]' '{print $3}' |awk -F'[_]' '{print $2}')
end=$(echo ${line} | awk -F'[#]' '{print $3}' |awk -F'[_]' '{print $3}')
fi
if [ -z "$end" ]; then
awk -v a=$seq -v b=$start -v c=$end -v d=$genome1 'BEGIN{OFS="\t"}  $1 == a {print d"#1#"$1,$2,$3,d"_"a"_"$2}' $anno
else
awk -v a=$seq -v b=$start -v c=$end -v d=$genome1 'BEGIN{OFS="\t"}  $1 == a && $2 >= b && $3 <= c {print d"#1#"$1"_"b"_"c,$2-b,$3-b,d"_"a"_"$2}' $anno
fi
done
