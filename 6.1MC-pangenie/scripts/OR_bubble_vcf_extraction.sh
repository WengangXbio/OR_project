cut -f1 combined.bubble.ORgt |sort |uniq | awk -F'_' 'BEGIN{OFS="\t"} {print $1"_"$2,$3,$3}' >  combined.bubble.ORgt.bed

for id in `ls *vcf |awk -F'.' '{print $1}'`; do
bedtools intersect -a ${id}.PanGenie_genotyping.vcf -b ../combined.bubble.ORgt.bed |cut -f1,2,10 |awk -F':' '{print $1}' > /share/home/yzwl_zhangwg/OR_project/MC-bos-31/pangenie/pangenie_output_or/${id}.PanGenie_genotyping.vcf.or
done



