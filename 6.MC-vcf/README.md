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
#All existing OR gene names (including reference-based and nonreference)
awk 'BEGIN{OFS="\t"} {sum=0; for(i=2; i<=32; i++) sum+=$i; print $0, sum}' ../or.PAV.matrix0.collection |awk '$33!=0 {print $1}' |sed '1d' |sed 's/^[^.]*\.//' > all_existing_OR.name
#Extracting all nonreference OR flanking sequences
grep -v "GCF_002263795" all_existing_OR.name | while IFS= read -r OR; do 
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
