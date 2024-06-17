chr=$1
start=$2
end=$3
OR=$4
coor=$(awk -v a=$OR '$4==a' all.func.bed)
genome=$(echo $coor |awk -F'#' '{print $1}')
cotig=$(echo $coor |awk -F'#' '{print $3}')
t1=$(echo $coor |awk '{print $2}')
t2=$(echo $coor |awk '{print $3}')
fasta=$(grep $genome genomelist |awk '{print $2}')
echo -e "$cotig\t$t1\t$t2\t$OR" |bedtools getfasta -name -fi ${fasta} -bed - > ${OR}_${start}_seqor.fa
forward=$(/public/home/qtyao/WZ_data_transfer/tools/seqtk-master/seqtk seq ${OR}_${start}_seqor.fa -U |tail -n1)
rc=$(/public/home/qtyao/WZ_data_transfer/tools/seqtk-master/seqtk seq ${OR}_${start}_seqor.fa -r -c -U|tail -n1)
awk -v a=$chr -v b=$start -v c=$end '$1==a && $2>b && $2<c' cattle31.raw.vcfbub.vcf > ${OR}_${start}_target.vcf
cln=$(($(grep -n $genome cattle31.raw.vcfbub.vcf.index |awk -F':' '{print $1}')+9))
cat ${OR}_${start}_target.vcf |while IFS= read -r vcf; do 
gt=$(echo $vcf |awk -v a=$cln '{print $a}')
alt=$(echo $vcf |awk '{print $5}')
ref=$(echo $vcf |awk '{print $4}')
vchr=$(echo $vcf |awk '{print $1}')
vpos=$(echo $vcf |awk '{print $2}')
if [[ $gt == 0 || $gt == "." ]]; then
echo -e "$vchr\t$vpos\t$ref"
else
valt=$(echo $alt |awk -F',' -v a=$gt '{print $a}')
echo -e "$vchr\t$vpos\t$valt"
fi
done > ${OR}_${start}_target.vcf.seq
uuf=$(grep $forward ${OR}_${start}_target.vcf.seq |cut -f1,2)
uur=$(grep $rc ${OR}_${start}_target.vcf.seq |cut -f1,2 )
if [ -n "$uuf" ]; then 
echo -e "$uuf\t$OR"
fi > ${OR}_${start}_nonrefOR_noproj.insertion_bubble
if [ -n "$uur" ]; then 
echo -e "$uur\t$OR"
fi >> ${OR}_${start}_nonrefOR_noproj.insertion_bubble
