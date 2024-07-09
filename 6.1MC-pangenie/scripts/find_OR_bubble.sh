cd /share/home/yzwl_zhangwg/OR_project/MC-bos-31/pangenie/openvcf
split -n l/500 ../cattle31.raw.vcfbub.vcf.diploid.noheader.removeN.vcf cattle31.raw.vcfbub.vcf.diploid.noheader.removeN.vcf

for suffix in `ls cattle31.raw.vcfbub.vcf.diploid.noheader.removeN.vcf* |awk -F'cattle31.raw.vcfbub.vcf.diploid.noheader.removeN.vcf' '{print $2}'|grep -v "^$"`; do
mkdir ${suffix}
cp cattle31.raw.vcfbub.vcf.diploid.noheader.removeN.vcf${suffix} all.coding.cdhit98.rename.fasta find_OR.sh ${suffix}/
cd ${suffix}
bsub -J findOR -n 1 -R "span[hosts=1]" -o %J.out -e %J.err -q cpu \
"
sh ./find_OR.sh ${suffix}
"
cd ..
done

find -name "*.ORgt"  -exec cat {} + > combined_ORgt.txt

