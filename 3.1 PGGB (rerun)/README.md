cd /share/home/yzwl_zhangwg/OR_project/PGGB/pan_cattle31_new
### Hybird MC 
```
for sequ in `seq 37328 37357`; do
mkdir NC_0${sequ}.1
for genome in `cat genomelist |grep -v "GCF_002263795"`; do
grep ${genome} MC-syntenic/cattle31.NC_0${sequ}.1.full.og.L |awk -v a=${genome} -F'[#:-]' 'BEGIN{OFS="\t"} {print $3,$5,$6,a"#"1"#"$3"_"$5"_"$6}'  |while IFS= read -r line; do
chr=$(echo ${line} |awk '{print $1}')
sta=$(echo ${line} |awk '{print $2}')
end=$(echo ${line} |awk '{print $3}')
name=$(echo ${line} |awk '{print $4}')
len=$(grep $chr /share/home/yzwl_zhangwg/OR_project/PGGB/assembly/${genome}.fasta.fai |awk '{print $2}')
if [ "$end" == "$len" ] && [ "$sta" == 0 ]; then
    echo -e "$chr\t$sta\t$end\t${genome}#1#$chr"
else
    echo -e "$chr\t$sta\t$end\t${name}"
fi
done |bedtools getfasta -fi /share/home/yzwl_zhangwg/OR_project/PGGB/assembly/${genome}.fasta -bed - -nameOnly > NC_0${sequ}.1/${genome}.fa
done 
samtools faidx /share/home/yzwl_zhangwg/OR_project/PGGB/assembly/GCF_002263795.fasta NC_0${sequ}.1 > NC_0${sequ}.1/GCF_002263795.fa
sed -i "s/NC_0${sequ}.1/GCF_002263795#1#NC_0${sequ}.1/g" NC_0${sequ}.1/GCF_002263795.fa
cat NC_0${sequ}.1/*.fa > NC_0${sequ}.1/NC_0${sequ}.1.fasta 
samtools faidx NC_0${sequ}.1/NC_0${sequ}.1.fasta 
chr=$((sequ - 37327))
cp NC_0${sequ}.1/NC_0${sequ}.1.fasta chr${chr}.pggb.fa
cp NC_0${sequ}.1/NC_0${sequ}.1.fasta.fai chr${chr}.pggb.fa.fai
done
```
