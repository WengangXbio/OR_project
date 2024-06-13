## MC vcf  
#### 1.

```
#Extract vcf site within OR regions 
bedtools intersect -wa -a cattle31.raw.vcfbub.vcf -b GCF_002263795.bed -u > cattle31.raw.vcfbub.vcf.OR
#output genotype of OR bubbles
cat cattle31.raw.vcfbub.vcf.OR | while IFS= read -r line; do 
chr=$(echo $line |awk '{print $1}')
str=$(echo $line |awk '{print $2}')
len=$(echo $line |awk '{print $4}' |awk '{print length($0)}')
end=$(($str+$len)); alt=$(echo $line |awk '{print $5}')
alt_len=$(echo $alt | awk -F, '{for(i=1;i<=NF;i++) printf "%d%s", length($i), (i<NF ? "," : "\n")}')
cor_or=$(echo -e "$chr\t$str\t$end" | bedtools intersect -b GCF_002263795.bed -a- -wb |awk '($3-$2)/($6-$5) ==1 {print $7}' |tr '\n' ',')
gt=$(echo $line |cut -d' ' -f10-)
if [ -n "$cor_or" ]; then 
echo -e "$chr\t$str\t$len\t$alt_len\t$cor_or\t$gt"
fi
done > bubbledOR.gt
#Bubbled OR names
cut -f5 bubbledOR.gt | sed 's/,/\n/g'  |awk 'NF' > bubbledOR.name
#Out-bubbled OR names
grep -v -f bubbledOR.name GCF_002263795.bed  > outbubbledOR.bed
#Out-bubbled OR genotype
cut -f4 outbubbledOR.bed | while IFS= read -r line; do
grep ${line} ../or_anno_transfer/GCF_002263795.func.bed > anno
chr=$(cat anno|awk -F "#" '{print $3}')
start=$(cat anno|awk '{print $2}')
end=$(cat anno|awk '{print $3}')
awk -v a=$chr -v b=$start -v c=$end '$1==a && $2>=b && $2<=c' cattle31.raw.vcfbub.vcf.OR |cut -f10- > genotyping
for ((j=1; j<=30; j++)); do
dot=$(awk -v col=$j '{print $col}' genotyping | grep '\.' |wc -l)
nodot=$(awk -v col=$j '{print $col}' genotyping | grep -v '\.' |wc -l)
echo -e "$dot","$nodot"
done | tr '\n' '\t' | awk -v a=${line} 'BEGIN{OFS="\t"} {print a,$0}'
done > outbubbledOR.gt
```
```
#All existing OR gene names (including reference-based and nonreference)
awk 'BEGIN{OFS="\t"} {sum=0; for(i=2; i<=32; i++) sum+=$i; print $0, sum}' ../or.PAV.matrix0.collection |awk '$33!=0 {print $1}' |sed '1d' |sed 's/^[^.]*\.//' > all_existing_OR.name
```
```
#All nonreference OR (including pseudo OR in reference or absent in reference)
grep -v "GCF_002263795" all_existing_OR.name |grep -f - or.PAV.matrix0.collection |awk '$28!=0  {print $1}' |sed 's/^[^.]*\.//' > ref_existing_pseudoOR.name
grep -v "GCF_002263795" all_existing_OR.name |grep -f - or.PAV.matrix0.collection |awk '$28==0  {print $1}' |sed 's/^[^.]*\.//' > ref_nonexisting_OR.name
```
```
#Extracting "ref_nonexisting_OR.name" flanking sequences
cat ref_nonexisting_OR.name | while IFS= read -r OR; do 
genome=$(cat ../all.func.bed |grep $OR |awk -F'#' '{print $1}')
chr=$(cat ../all.func.bed |grep $OR |awk -F'#' '{print $3}')
start=$(($(cat ../all.func.bed |grep $OR |awk '{print $2}') - 50000))
end=$(($(cat ../all.func.bed |grep $OR |awk '{print $3}')+ 50000))
if [ "$start" -lt 0 ]; then 
start=0
fi
fasta=$(grep $genome ../../genomelist |awk '{print $2}')
echo -e "$chr\t$start\t$end\t$OR" |bedtools getfasta -name -fi ${fasta} -bed - 
done > all_existing_OR.nonref.fasta
fasta=/public/home/qtyao/WZ_data_transfer/cactus/bos/assembly/GCF_002263795.3_ARS-UCD2.0_genomic.fna
minimap2 -cx asm5 ${fasta} all_existing_OR.nonref.fasta  -t 4 >  all_existing_OR.nonref.paf
```
```
#project nonreference OR on reference
for OR in `cat ref_nonexisting_OR.name ` ; do
chr=$(grep ${OR} or.PAV.matrix0.collection |awk -F'.' '{print $1}')
cotig=$(awk -v a=${chr} '$1==a' ../seq.txt |cut -f2)
start=$(($(awk -v a=$OR -v b=$cotig '$1==a && $6==b' all_existing_OR.nonref.paf |cut -f8,9 |tr '\t' '\n' |sort -k1,1n |head -n1)-10000))
end=$(($(awk -v a=$OR -v b=$cotig '$1==a && $6==b' all_existing_OR.nonref.paf |cut -f8,9 |tr '\t' '\n' |sort -k1,1n |tail -n1)+10000))
echo -e "$cotig\t$start\t$end\t$OR" 
done > ref_nonexisting_OR.name.insertion.bed
awk '$2!=-10000 && $3-$2<200000'  ref_nonexisting_OR.name.insertion.bed >  ref_nonexisting_OR.name.insertion.bed.filter
```
```
cat all_existing_OR.nonref.insertion.bed.filter | while IFS= read -r line; do 
chr=$(echo $line |awk '{print $1}')
start=$(echo $line |awk '{print $2}')
end=$(echo $line |awk '{print $3}')
OR=$(echo $line |awk '{print $4}')
coor=$(awk -v a=$OR '$4==a' all.func.bed)
genome=$(echo $coor |awk -F'#' '{print $1}')
cotig=$(echo $coor |awk -F'#' '{print $3}')
t1=$(echo $coor |awk '{print $2}')
t2=$(echo $coor |awk '{print $3}')
fasta=$(grep $genome genomelist |awk '{print $2}')
echo -e "$cotig\t$t1\t$t2\t$OR" |bedtools getfasta -name -fi ${fasta} -bed - > seqor.fa
forward=$(/public/home/qtyao/WZ_data_transfer/tools/seqtk-master/seqtk seq seqor.fa -U |tail -n1)
rc=$(/public/home/qtyao/WZ_data_transfer/tools/seqtk-master/seqtk seq seqor.fa -r -c -U|tail -n1)
awk -v a=$chr -v b=$start -v c=$end '$1==a && $2>b && $2<c' cattle31.raw.vcfbub.vcf > target.vcf
cln=$(($(grep -n $genome cattle31.raw.vcfbub.vcf.index |awk -F':' '{print $1}')+9))
cat target.vcf |while IFS= read -r vcf; do 
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
done > target.vcf.seq
uuf=$(grep $forward target.vcf.seq |cut -f1,2)
uur=$(grep $rc target.vcf.seq |cut -f1,2 )
if [ -n "$uuf" ]; then 
echo -e "$uuf\t$OR"
fi
if [ -n "$uur" ]; then 
echo -e "$uur\t$OR"
fi
done 

```


