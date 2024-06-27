input_vcf=$1
input_node=$2
reference_id=$3
grep "#CHROM" ${input_vcf} |cut -f10- |tr '\t' '\n' > vcf.index
grep -v '^#' ${input_vcf} |while IFS= read -r line; do
CHR=$(echo $line |awk '{print $1}')
POS=$(echo $line |awk '{print $2}')
VS=$(echo $line |awk '{print $8}' |awk -F 'VS=' '{print $2}' |awk -F';' '{print $1}' |tr -d '>')
AWALK=$(echo $line |awk '{print $8}' |awk -F 'AWALK=' '{print $2}')
poly=$(echo $AWALK |tr ',' '\n' |wc -l)
if [ $poly -eq 1 ]; then
continue
fi
gt=$(echo $line |cut -d' ' -f10-)
tg=$(echo $gt |tr ' ' '\n'| awk -F':' '{print $1}')
VSlast=$(awk -v a=$VS '$2==a {print $3}'  ${input_node} |rev |cut -c 1)
alt_seq=$(echo "${AWALK}" |tr ',' '\n' |while IFS= read -r alt; do
echo "$alt" | sed 's/>\|</\n& /g' |grep -v "^$" |while IFS= read -r node; do
strand=$(echo "$node" |awk '{print $1}')
sss=$(echo "$node" |awk '{print $2}')
if [[ ${strand} == ">" ]];then
awk -v a=$sss '$2==a {print $3}'  ${input_node}
fi
if [[ ${strand} == "<" ]];then
awk -v a=$sss '$2==a {print $3}'  ${input_node} |awk '{print ">read\n"$0}' |~/WZ_data_transfer/tools/seqtk-master/seqtk seq -r |sed '1d'
fi
if [[ ${strand} == "*" ]];then
echo "ZW"
fi
done |tr -d '\n' |awk '{print $0","}'
done |tr -d '\n' |sed 's/,$//'|tr ',' '\n'|awk -v a=${VSlast} '{print a$0}' |tr '\n' ',' |tr -d 'ZW' |sed 's/,$//')
echo $alt_seq |tr ',' '\n' | awk '{print NR-1,$0}' > temp1
echo "." >> temp1
cat temp1 | while IFS= read -r reld; do
ttgg=$(echo $reld |awk -F' ' '{print $1}')
id=$(echo $tg |tr ' ' '\n' |paste vcf.index - |awk -v a=${ttgg} '$2==a {print $1}'  |tr '\n' ',' )
echo -e "$reld $id"
done > temp2
grep ${reference_id} temp2 > temp3
grep -v ${reference_id} temp2 >> temp3
awk '{ if ($1 != ".") {print NR-1,$2,$3}else{print $0}}' temp3 > temp4
outGT=$(cat vcf.index | while IFS= read -r ind; do
grep $ind temp4 |awk '{print $1}'
done |tr '\n' '\t')
REF=$(head -1 temp4 |cut -d' ' -f2)
ALT=$(sed '1d' temp4 |awk '$1 != "."'|cut -d' ' -f2 |tr '\n' ','|sed 's/,$//')
echo -e "$CHR\t$POS\t$CHR"_"$POS\t$REF\t$ALT\t60\t"."\tAC=1;AF=0;AN=4;NS=4\tGT\t$outGT"
done > vcf.body
