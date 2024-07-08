vcfbub -l 0 -r 100000 --input cattle31.raw.vcf  > cattle31.raw.vcfbub.vcf 

grep -v "^#" cattle31.raw.vcfbub.vcf |cut -f10- > cattle31.raw.vcfbub.vcf.gt
sed -r 's/([0-9]+|\.)/&|\1/g'  cattle31.raw.vcfbub.vcf.gt > cattle31.raw.vcfbub.vcf.gt.diploid
grep -v "^#" cattle31.raw.vcfbub.vcf |cut -f1-9 |awk -F '\t' 'BEGIN{OFS="\t"} {print $1,$2,".",$4,$5,".","PASS",".","GT"}' > cattle31.raw.vcfbub.vcf.site_info
paste -d '\t' cattle31.raw.vcfbub.vcf.site_info cattle31.raw.vcfbub.vcf.gt.diploid > cattle31.raw.vcfbub.vcf.bottom
grep "^#" cattle31.raw.vcfbub.vcf > cattle31.raw.vcfbub.vcf.top
cat  cattle31.raw.vcfbub.vcf.top cattle31.raw.vcfbub.vcf.bottom > cattle31.raw.vcfbub.vcf.diploid
