### 1.find OR's bubble in VCF file
```
split -n l/999 cattle31.raw.vcfbub.vcf cattle31.raw.vcfbub.vcf

for suffix in `ls cattle31.raw.vcfbub.vcf* |awk -F'cattle31.raw.vcfbub.vcf' '{print $2}'|grep -v "^$"`; do
mkdir ${suffix}
cp cattle31.raw.vcfbub.vcf${suffix} all.coding.cdhit98.rename.fasta find_OR.sh ${suffix}/
cd ${suffix}
bsub -J findOR -n 2 -R "span[hosts=1]" -o %J.out -e %J.err -q cpu \
"
sh ./find_OR.sh ${suffix}
"
cd ..
done
```

### 2. filter VCF with OR SV
```
cut -f1 combined.bubble.ORgt |sort |uniq | awk -F'_' 'BEGIN{OFS="\t"} {print $1"_"$2,$3,$3}' >  combined.bubble.ORgt.bed
for id in `ls *vcf |awk -F'.' '{print $1}'`; do
bedtools intersect -a ${id}.PanGenie_genotyping.vcf -b /share/home/yzwl_zhangwg/OR_project/MC-bos-31/openvcf/combined.bubble.ORgt.bed |cut -f1,2,10 |awk -F':' '{print $1}' >/share/home/yzwl_zhangwg/OR_project/MC-bos-31/filter_output/${id}.PanGenie_genotyping.vcf.or
done
```
