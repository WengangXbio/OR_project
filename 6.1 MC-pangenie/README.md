## find OR's bubble in VCF file
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

