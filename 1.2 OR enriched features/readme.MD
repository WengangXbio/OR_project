## Genome-wide profiling based on pangenome and featured with OR regions

### window depth by 10K window project on ARS2.0 
```
for nchr in `seq 1 29`  ; do
odgi build -t 4 -P -g chr${nchr}.gfa -o chr${nchr}.og -O
rhe=$(echo chr${nchr}) 
che=$(awk -v a=${rhe} '$1==a {print $2}' chr_anno.txt)
odgi depth -i ${rhe}.og -r ${che}#1#${rhe} | bedtools makewindows -b /dev/stdin -w 10000 > window_depth/${rhe}.w10k.bed
odgi depth -i ${rhe}.og -b window_depth/${rhe}.w10k.bed --threads 2 | bedtools sort > window_depth/${rhe}.w10k.depth.bed
done
```

### variants density on CDS regions by 10K window walking
```
awk -F'\t' 'BEGIN{OFS="\t"} $3=="CDS" {print $1,$4,$5}' GCF_002263795.3_ARS-UCD2.0_genomic.gff > GCF_002263795.3_ARS-UCD2.0_genomic.gff.CDS.bed
for nchr in `seq 1 29`  ; do
chr="chr${nchr}"
ass=$(awk -v a=${chr} '$1==a {print $2}' chr_anno.txt)
sed -i "s/${ass}/${chr}/g" GCF_002263795.3_ARS-UCD2.0_genomic.gff.CDS.bed
done
grep "chr" GCF_002263795.3_ARS-UCD2.0_genomic.gff.CDS.bed > GCF_002263795.3_ARS-UCD2.0_genomic.gff.CDS.bed0
for nchr in `seq 1 29`  ; do
chr="chr${nchr}"
length=$(awk -v a=${chr} '$1==a {print $3}' chr_anno.txt)
grep -v "^#" /share/home/yzwl_zhangwg/OR_project/PGGB/work/${chr}.vcf | awk -v a=${chr} 'BEGIN{OFS="\t"} {print a,$2,$2+1}' - > tmp.vcf.bed
bedtools intersect -b tmp.vcf.bed -a GCF_002263795.3_ARS-UCD2.0_genomic.gff.CDS.bed0 > tmp.CDS.vcf.bed
echo -e "$chr\t"0"\t$length" |bedtools makewindows -b - -w 1000 > tmp.window.bed
bedtools intersect -b tmp.CDS.vcf.bed -a tmp.window.bed -c > window_vcf_density/${chr}.CDS.vcf.density
rm tmp.vcf.bed tmp.CDS.vcf.bed tmp.window.bed
done
```

### variants density by 10K window walking
```
for nchr in `seq 1 29`  ; do
chr="chr${nchr}"
length=$(awk -v a=${chr} '$1==a {print $3}' chr_anno.txt)
grep -v "^#" /share/home/yzwl_zhangwg/OR_project/PGGB/work/${chr}.vcf | awk -v a=${chr} 'BEGIN{OFS="\t"} {print a,$2,$2+1}' - > tmp.vcf.bed
echo -e "$chr\t"0"\t$length" |bedtools makewindows -b - -w 10000 > tmp.window.bed
bedtools intersect -b tmp.vcf.bed -a tmp.window.bed -c > window_vcf_density/${chr}.vcf.w10k.density
rm tmp.window.bed tmp.vcf.bed
done
```

### PAV average value by 10K window walking
```
cp /share/home/yzwl_zhangwg/OR_project/PGGB/work/window_depth/*.w10k.bed ./window_pav
cd /share/home/yzwl_zhangwg/OR_project/PGGB/work/
for nchr in `seq 1 29`  ; do
chr="chr${nchr}"
odgi pav -i ${chr}.og -b window_pav/${chr}.w10k.bed > window_pav/${chr}.w10k.pav
done
```
