cd /share/home/yzwl_zhangwg/OR_project/MC-bos-31/pangenie/removeN
grep -v "^#" cattle31.raw.vcfbub.vcf.diploid > cattle31.raw.vcfbub.vcf.diploid.noheader
split -n l/150 ../cattle31.raw.vcfbub.vcf.diploid.noheader cattle31.raw.vcfbub.vcf.diploid.noheader
for vcf in `cat list `; do
bsub -J vcf -n 1 -R "span[hosts=1]" -o %J.out -e %J.err -q cpu \
"
./removeN.sh ${vcf}
"
done
cat *.removeN > cattle31.raw.vcfbub.vcf.diploid.noheader.removeN
grep "^#" ../cattle31.raw.vcfbub.vcf.diploid |cat - cattle31.raw.vcfbub.vcf.diploid.noheader.removeN > ../cattle31.raw.vcfbub.vcf.diploid.noheader.removeN.vcf

#cat removeN.sh 
input=$1
mkdir ${input}_temp
echo "" > ${input}.NNN.txt
cat ${input} | while IFS= read -r line; do
Ngt=$(echo $line |awk '{print $5}' |tr ',' '\n' |grep -n "N" |awk -F':' '{print $1}')
if [ -n "$Ngt" ]; then
gt=$(echo $line |cut -d ' ' -f10-)
echo $line |awk '{print $5}' |tr ',' '\n' |awk '{print NR,NR,$0}' > ${input}_temp/alt_seq
for agt in `echo $Ngt` ; do
gt=$(echo $gt |tr ' ' '\n' |awk -F'|' -v a=$agt 'BEGIN{OFS=""}{if ($1==a) $1="."; if ($2==a) $2="."; print $1,"|",$2}')
awk -v a=$agt '$1!=a' ${input}_temp/alt_seq > ${input}_temp/temp
mv ${input}_temp/temp ${input}_temp/alt_seq
done
echo $line |awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' > ${input}_temp/c1234_txt
awk '{print $3}' ${input}_temp/alt_seq |tr '\n' ',' |sed 's/,$//' |awk '{print $0'\n'}' > ${input}_temp/c5_txt
echo $line |awk 'BEGIN{OFS="\t"} {print $6,$7,$8,$9}' > ${input}_temp/c6789_txt
echo $gt |tr ' ' '\t' > ${input}_temp/gt_txt
awk '$2=NR {print $1,$2}' ${input}_temp/alt_seq |while IFS= read -r yty; do
old=$(echo $yty |awk '{print $1}')
new=$(echo $yty |awk '{print $2}')
cat ${input}_temp/gt_txt | tr '\t' '\n' | awk -F'|' -v a="$old" -v b="$new" 'BEGIN{OFS=""}{if ($1==a) $1=b; if ($2==a) $2=b; print $1,"|",$2}' > ${input}_temp/temp_gt
cat ${input}_temp/temp_gt | tr '\n' '\t' > ${input}_temp/gt_txt
rm ${input}_temp/temp_gt
done
lnl=$(cat ${input}_temp/c5_txt |wc -l )
if [ "$lnl" -gt 0 ]; then
paste -d '\t' ${input}_temp/c1234_txt ${input}_temp/c5_txt ${input}_temp/c6789_txt ${input}_temp/gt_txt
echo  $line |awk '{print $1,$2}' >> ${input}.NNN.txt
else
echo  $line |awk '{print $1,$2,"revomal"}' >> ${input}.NNN.txt
fi
rm ${input}_temp/c1234_txt ${input}_temp/c5_txt ${input}_temp/c6789_txt ${input}_temp/gt_txt ${input}_temp/alt_seq
else
echo $line |sed 's/ /\t/g'
fi
done > ${input}.removeN
