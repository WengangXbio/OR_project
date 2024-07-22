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

### Run PGGB
```
cd /share/home/yzwl_zhangwg/OR_project/PGGB/pan_cattle31_new
sh pggb.sh
```

### 1. Creating OR annotation by adjusting coordinates
```
cd /share/home/yzwl_zhangwg/OR_project/PGGB/pan_cattle31_new/work
for nchr in `seq 1 30` ; do
CHROM=chr${nchr}
mkdir $CHROM && cd $CHROM 
ln -s /share/home/yzwl_zhangwg/OR_project/PGGB/pan_cattle31_new/scripts/coding_annotation_generator.sh ./
ln -s /share/home/yzwl_zhangwg/OR_project/PGGB/pan_cattle31_new/scripts/dead_annotation_generator.sh ./
ln -s /share/home/yzwl_zhangwg/OR_project/PGGB/pan_cattle31_new/scripts/precess_bed.sh ./
ln -s /share/home/yzwl_zhangwg/OR_project/PGGB/pan_cattle31_new/scripts/modify_coor.sh ./
chmod +x modify_coor.sh
chmod +x precess_bed.sh
chmod +x coding_annotation_generator.sh
chmod +x dead_annotation_generator.sh
fastapath=/share/home/yzwl_zhangwg/OR_project/PGGB/pan_cattle31_new
OR_path=/share/home/yzwl_zhangwg/OR_project/genome2OR
./coding_annotation_generator.sh ${fastapath} ${OR_path} ${CHROM}
./dead_annotation_generator.sh ${fastapath} ${OR_path} ${CHROM}
sed -i "s/OMG/${og}/g" projection.sh
csub < projection.sh
cd ..
done
```

### 2. Adjusting OR projections by collapsing coordinates
```
for nchr in `seq 1 30` ; do
cd chr${nchr}
if [ -s "./injref.untangle.tsv" ]; then
cp /share/home/yzwl_zhangwg/OR_project/PGGB/pan_cattle31_new/scripts/count_run.sh ./
sed -i "s/CHRN/${nchr}/g" count_run.sh
chmod +x count_run.sh
csub < count_run.sh
else 
   echo "chr${nchr} is blank"
fi
cd ..
done
```

### 2.1 Collapsing coordinates is slow when OR is highly enriched in certain chromosome (e.g. chromosome 7 and 15), running large collasping scripts as follow
```
#1.run first
for iseq in `cat inputseq` ; do
mkdir ${iseq}_collaspe 
cd ${iseq}_collaspe 
cp /share/home/yzwl_zhangwg/OR_project/PGGB/pan_cattle31_new/scripts/collaspe_untangle.sh collaspe_untangle.sh 
cp /share/home/yzwl_zhangwg/OR_project/PGGB/pan_cattle31_new/scripts/large_collaspe.sh large_collaspe.sh
cp ../injref.untangle.tsv ./
echo $iseq > inputseq
csub < large_collaspe.sh
cd ..
done

#2.run first
echo > collection 
for splitfile in `ls -l ./ | grep '^d' | grep "_collaspe" | awk '{print $9}'`; do
cat collection ${splitfile}/injref.untangle.tsv.collapse > temp
mv temp collection
done
mv collection injref.untangle.tsv.collapse
```


