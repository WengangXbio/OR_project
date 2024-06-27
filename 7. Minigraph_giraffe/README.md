## Minigraph + giraffe 
Using minigraph build graph and then mapping short reads to call SV
### Seperate sequence by chromosome
```
cd /public/home/qtyao/WZ_data_transfer/minigraph/1.minigraph
for ctg in `seq 37328 37357 |awk '{print} END {print 82638}'` ; do
mkdir NC_0${ctg}.1
cd  NC_0${ctg}.1
for genome in `ls /public/home/qtyao/WZ_data_transfer/minigraph/0.assembly/ |grep fasta |grep -v fai |awk -F'.' '{print $1}'` ; do
if [[ ${genome} != "GCF_002263795" ]]; then
awk -F'[#:-]' -v a=${genome} 'BEGIN{OFS="\t"} $1==a {print $3,$5,$6,$3"_"$5}' ../cattle31.NC_0${ctg}.1.syntenic.regions | \
bedtools getfasta -fi /public/home/qtyao/WZ_data_transfer/minigraph/0.assembly/${genome}.fasta  -bed - -nameOnly > ${genome}_NC_0${ctg}.1.fasta
fi
if [[ ${genome} == "GCF_002263795" ]]; then
awk -F'[#:-]' -v a=${genome} 'BEGIN{OFS="\t"} $1==a {print $3,$4,$5,$3"_"$4}' ../cattle31.NC_0${ctg}.1.syntenic.regions | \
bedtools getfasta -fi /public/home/qtyao/WZ_data_transfer/minigraph/0.assembly/${genome}.fasta  -bed - -nameOnly > ${genome}_NC_0${ctg}.1.fasta
fi
done
cd /public/home/qtyao/WZ_data_transfer/minigraph/1.minigraph
done
```
### Implement minigraph by chromosome
```
for ctg in `seq 37328 37357 |awk '{print} END {print 82638}'` ; do
cd  NC_0${ctg}.1
chr=NC_0${ctg}.1
awk -v a=${chr} '{print $0"_"a".fasta"}' ../genome_prefix |tr '\n' ' ' |\
awk -v a=$chr '{print "/public/home/qtyao/WZ_data_transfer/tools/minigraph-master/minigraph -cxggs -t24 "$0" > "a".gfa" }' > minigraph.sh
bsub -J minigraph -n 24 -R "span[hosts=1] rusage[mem=6GB] " -o %J.out -e %J.err -q normal \
'
sh minigraph.sh
'
cd /public/home/qtyao/WZ_data_transfer/minigraph/1.minigraph
done
```

### Call SV on each assembly
```
for ctg in `seq 37328 37357 ` ; do
cd  NC_0${ctg}.1
chr=NC_0${ctg}.1
gfa=${chr}.gfa
for genome in `ls /public/home/qtyao/WZ_data_transfer/minigraph/0.assembly/ |grep fasta |grep -v fai |awk -F'.' '{print $1}'` ; do
bsub -J bed -n 12 -R "span[hosts=1] rusage[mem=2GB] " -o %J.out -e %J.err -q normal \
"
/public/home/qtyao/WZ_data_transfer/tools/minigraph-master/minigraph -cxasm --call -t12 ${gfa} ${genome}_NC_0${ctg}.1.fasta > ${genome}.bed
"
done
cd /public/home/qtyao/WZ_data_transfer/minigraph/1.minigraph
done
```
### Merge SV for all assemblies to generate minigraph-based VCF
```
for ctg in `seq 37328 37357 ` ; do
cd  NC_0${ctg}.1
chr=NC_0${ctg}.1
awk '{print $0".bed"}' ../genome_prefix |tr '\n' ' ' |xargs paste | /public/home/qtyao/WZ_data_transfer/tools/mg-cookbook-v1_x64-linux/k8 \
/public/home/qtyao/WZ_data_transfer/tools/mg-cookbook-v1_x64-linux/mgutils.js merge -s <(cat ../genome_prefix ) - | gzip > NC_0${ctg}.1.sv.bed.gz
/public/home/qtyao/WZ_data_transfer/tools/mg-cookbook-v1_x64-linux/k8 /public/home/qtyao/WZ_data_transfer/tools/mg-cookbook-v1_x64-linux/mgutils-es6.js \
merge2vcf -r0 NC_0${ctg}.1.sv.bed.gz > NC_0${ctg}.1.sv.vcf
grep "^##" NC_0${ctg}.1.sv.vcf | cat - ../replace_title - <(grep -v "^#" NC_0${ctg}.1.sv.vcf) > NC_0${ctg}.1.sv.vcf0
cd /public/home/qtyao/WZ_data_transfer/minigraph/1.minigraph
done
```

### Modifing vcf for each chromosome
```
for ctg in `seq 37328 37357 ` ; do
chr=NC_0${ctg}.1
mkdir ${chr}
cd ${chr}
ln -s /public/home/qtyao/WZ_data_transfer/minigraph/1.minigraph/${chr}/${chr}.sv.vcf0 ${chr}.sv.vcf0
ln -s /public/home/qtyao/WZ_data_transfer/minigraph/1.minigraph/${chr}/${chr}.gfa ${chr}.gfa 
ln -s /public/home/qtyao/WZ_data_transfer/minigraph/2.vgvcf/vcf_generation.sh vcf_generation.sh
grep '^S' ${chr}.gfa > out.gfa.node
bsub -J bed -n 2 -R "span[hosts=1] rusage[mem=4GB] " -o %J.out -e %J.err -q normal \
"
sh ./vcf_generation.sh ${chr}.sv.vcf0 out.gfa.node GCF_002263795  # input_vcf, input_node, reference_id
"
cd /public/home/qtyao/WZ_data_transfer/minigraph/2.vgvcf
done
```

### Correcting chromosome name of modified vcf
```
for ctg in `seq 37328 37357 ` ; do
chr=NC_0${ctg}.1
chr0=NC_0${ctg}.1_0
cd ${chr}
sed "s/${chr0}/${chr}/g" vcf.body >  ${chr}.sv.vcf0.vg.vcf0
cut -f1 ${chr}.sv.vcf0.vg.vcf0 |uniq
cd /public/home/qtyao/WZ_data_transfer/minigraph/2.vgvcf
done
```

### Find inconsistent site (including inconsistent sites and duplicated nodes)
_Inconsistent sites_: The modified vcf's reference seqence is not consistent with reference sequence

_Duplicated nodes_: Some of sequences of two alternative genotypes are exactly same  
```
for ctg in `seq 37328 37357` ; do
chr=NC_0${ctg}.1
cd ${chr}
echo > conflict_genotype
echo > inconsistent_with_ref
cp ../vcf_correction.sh ./
bsub -J correction -n 2 -R "span[hosts=1] rusage[mem=4GB] " -o %J.out -e %J.err -q normal \
"
sh ./vcf_correction.sh ${chr}
"
cd /public/home/qtyao/WZ_data_transfer/minigraph/2.vgvcf
done
```

### Merge vcf
```
cat header \
./NC_037328.1/NC_037328.1.sv.vcf0.vg.correction.vcf0 \
./NC_037329.1/NC_037329.1.sv.vcf0.vg.correction.vcf0 \
./NC_037330.1/NC_037330.1.sv.vcf0.vg.correction.vcf0 \
./NC_037331.1/NC_037331.1.sv.vcf0.vg.correction.vcf0 \
./NC_037332.1/NC_037332.1.sv.vcf0.vg.correction.vcf0 \
./NC_037333.1/NC_037333.1.sv.vcf0.vg.correction.vcf0 \
./NC_037334.1/NC_037334.1.sv.vcf0.vg.correction.vcf0 \
./NC_037335.1/NC_037335.1.sv.vcf0.vg.correction.vcf0 \
./NC_037336.1/NC_037336.1.sv.vcf0.vg.correction.vcf0 \
./NC_037337.1/NC_037337.1.sv.vcf0.vg.correction.vcf0 \
./NC_037338.1/NC_037338.1.sv.vcf0.vg.correction.vcf0 \
./NC_037339.1/NC_037339.1.sv.vcf0.vg.correction.vcf0 \
./NC_037340.1/NC_037340.1.sv.vcf0.vg.correction.vcf0 \
./NC_037341.1/NC_037341.1.sv.vcf0.vg.correction.vcf0 \
./NC_037342.1/NC_037342.1.sv.vcf0.vg.correction.vcf0 \
./NC_037343.1/NC_037343.1.sv.vcf0.vg.correction.vcf0 \
./NC_037344.1/NC_037344.1.sv.vcf0.vg.correction.vcf0 \
./NC_037345.1/NC_037345.1.sv.vcf0.vg.correction.vcf0 \
./NC_037346.1/NC_037346.1.sv.vcf0.vg.correction.vcf0 \
./NC_037347.1/NC_037347.1.sv.vcf0.vg.correction.vcf0 \
./NC_037348.1/NC_037348.1.sv.vcf0.vg.correction.vcf0 \
./NC_037349.1/NC_037349.1.sv.vcf0.vg.correction.vcf0 \
./NC_037350.1/NC_037350.1.sv.vcf0.vg.correction.vcf0 \
./NC_037351.1/NC_037351.1.sv.vcf0.vg.correction.vcf0 \
./NC_037352.1/NC_037352.1.sv.vcf0.vg.correction.vcf0 \
./NC_037353.1/NC_037353.1.sv.vcf0.vg.correction.vcf0 \
./NC_037354.1/NC_037354.1.sv.vcf0.vg.correction.vcf0 \
./NC_037355.1/NC_037355.1.sv.vcf0.vg.correction.vcf0 \
./NC_037356.1/NC_037356.1.sv.vcf0.vg.correction.vcf0 \
./NC_037357.1/NC_037357.1.sv.vcf0.vg.correction.vcf0 > all.sv.vcf0.vg.correction.vcf0
```

### Preparation of zipped vcf for each chromosome
```
/home/zhangwengang/tools/bcftools-1.2/htslib-1.2.1/bgzip -c all.sv.vcf0.vg.correction.vcf0 > all.sv.vcf0.vg.correction.vcf0.gz
/home/zhangwengang/tools/bcftools-1.2/htslib-1.2.1/tabix -p vcf all.sv.vcf0.vg.correction.vcf0.gz
for ctg in `seq 37328 37357` ; do
chr=NC_0${ctg}
grep "^#\|^${chr}" all.sv.vcf0.vg.correction.vcf0 > chunk_${chr}.1.correction.vcf0
/home/zhangwengang/tools/bcftools-1.2/htslib-1.2.1/bgzip -c chunk_${chr}.1.correction.vcf0 > chunk_${chr}.1.correction.vcf0.gz
/home/zhangwengang/tools/bcftools-1.2/htslib-1.2.1/tabix -p vcf chunk_${chr}.1.correction.vcf0.gz
done
```

### Index graph
```
vg construct -f -a -S -r GCF_002263795.fasta -v  all.sv.vcf0.vg.correction.vcf0.gz -p > out.vg
vg index out.vg -x out.xg  -t 30 -L -p --temp-dir tmp1
vg gbwt -x out.vg -v all.sv.vcf0.vg.correction.vcf0.gz -o out.gbwt -p --temp-dir tmp1
vg gbwt -x out.vg -E -o out.paths.gbwt -p --temp-dir tmp1
vg gbwt -m out.gbwt out.paths.gbwt -o out.combined.gbwt -p --temp-dir tmp1
vg gbwt -x out.vg out.combined.gbwt --gbz-format -g out.gbz -p --temp-dir tmp1
vg index out.gbz -j out.dist  -t 30 --temp-dir tmp1
vg minimizer out.gbz -d out.dist -o out.min
vg snarls out.xg > out.snarls 
```
### Giraffle mapping
```
vg giraffe -p -t 36 -Z out.gbz -d out.dist -m out.min \
    -f SRR12094761_1.clean.fq.gz -f SRR12094761_2.clean.fq.gz > SRR12094761.gam
vg chunk -t 30 -x out.vg -a SRR12094761.gam -M -g 

#calling
#vg pack -x out.vg -g chunk_NC_037328.1.gam -Q 5 -t 48 -o chunk_NC_037328.1.pack 
#vg call out.vg -k NC_037328.1.pack --vcf chunk_NC_037328.1.correction.vcf0.gz -r out.snarls -t 48 > NC_037328.1.vcf

seq 37348 37357 | parallel -j 14 '
  vg pack -x out.vg -g chunk_NC_0{}.1.gam -Q 5 -t 14 -o chunk_NC_0{}.1.pack &&
  vg call out.vg -k chunk_NC_0{}.1.pack --vcf chunk_NC_0{}.1.correction.vcf0.gz -r out.snarls -t 14 > chunk_NC_0{}.1.vcf
'
```

### Preparing funcational OR annotation for all assemblies (suitable for adjusted coordinates when trim seqences for chromosomal graph construction)
```
for genome in `ls /public/home/qtyao/WZ_data_transfer/minigraph/0.assembly/ |grep fasta |grep -v fai |awk -F'.' '{print $1}'` ; do
if [[ ${genome} != "GCF_002263795" ]]; then
anno_path=/public/home/qtyao/WZ_data_transfer/cactus/bos/work/or_anno_transfer
coding_anno=$(ls ${anno_path} |grep ${genome} |grep g2or_itera2_func.bed |grep -v dead)
cat /public/home/qtyao/WZ_data_transfer/minigraph/1.minigraph/*.syntenic.regions | awk -F'[#:-]' -v a=${genome} 'BEGIN{OFS="\t"} $1==a {print $3,$5,$6,$3"_"$5}' |\
bedtools intersect -a - -b ${anno_path}/${coding_anno} -wa -wb |\
awk -v a=${genome} 'BEGIN{OFS="\t"} {print $4,$2+$6,$2+$7,a"_"$5"_"$6}' > ${genome}.func.bed
fi
done
genome=GCF_002263795
coding_anno=$(ls ${anno_path} |grep ${genome} |grep g2or_itera2_func.bed |grep -v dead)
awk -v a=${genome} 'BEGIN{OFS="\t"} {print $1,$2,$3,a"_"$1"_"$2}' ${anno_path}/${coding_anno} >${genome}.func.bed
```

### Extracting reference bubble OR VCF and summarizing (ref_bubble.OR.vcf.gt)
```
bedtools intersect -wa -a all.sv.vcf0.vg.correction.vcf0 -b anno/GCF_002263795.func.bed -u > ref_bubble.OR.vcf
cat ref_bubble.OR.vcf | while IFS= read -r line; do 
chr=$(echo $line |awk '{print $1}')
str=$(echo $line |awk '{print $2}')
len=$(echo $line |awk '{print $4}' |awk '{print length($0)}')
end=$(($str+$len)); alt=$(echo $line |awk '{print $5}')
alt_len=$(echo $alt | awk -F, '{for(i=1;i<=NF;i++) printf "%d%s", length($i), (i<NF ? "," : "\n")}')
cor_or=$(echo -e "$chr\t$str\t$end" | bedtools intersect -b GCF_002263795.bed -a - -wb |awk '($3-$2)/($6-$5) ==1 {print $7}' |tr '\n' ',')
gt=$(echo $line |cut -d' ' -f10-)
if [ -n "$cor_or" ]; then 
echo -e "$chr\t$str\t$len\t$alt_len\t$cor_or\t$gt"
fi
done > ref_bubble.OR.vcf.gt
```

### Localization of non-reference OR's coordinates corresponding to reference vcf (intermedian file)
```
for col in `seq 7 36` ; do
line=$((col -5))
genome=$(awk -v a=$line 'NR==a' ../1.minigraph/genome_prefix)
for ctg in `seq 37328 37357` ; do
chr=NC_0${ctg}.1
zcat ../1.minigraph/${chr}/${chr}.sv.bed.gz | grep -v "^#" |\
awk -v a=$chr -v b=${col} '{print a"_"$2":"$b}' |awk -F':' 'BEGIN{OFS="\t"} $4!="" {print $4,$5,$6,$1}' 
done > nonref_bubble_coordinate/${genome}.bed
done 

anno_path=/public/home/qtyao/WZ_data_transfer/cactus/bos/work/or_anno_transfer
for genome in `ls /public/home/qtyao/WZ_data_transfer/minigraph/0.assembly/ |grep fasta |grep -v fai |awk -F'.' '{print $1}' |grep -v GCF_002263795` ; do
echo $genome
sh ./modify_coor.sh nonref_bubble_coordinate/${genome}.bed nonref_bubble_coordinate/${genome}.bed0
bedtools sort -i nonref_bubble_coordinate/${genome}.bed0 > nonref_bubble_coordinate/${genome}.sort.bed0
bedtools merge -i nonref_bubble_coordinate/${genome}.sort.bed0 -c 4 -o collapse > nonref_bubble_coordinate/${genome}.sort.collapse.bed0
bedtools intersect -a anno/${genome}.func.bed -b nonref_bubble_coordinate/${genome}.sort.collapse.bed0 -wa -wb > nonref_bubble_coordinate/${genome}.sort.collapse.interectOR.bed0 
awk '{print $1"_"$2}' ref_bubble.OR.vcf.gt |grep -f - nonref_bubble_coordinate/${genome}.sort.collapse.interectOR.bed0 |\
cut -f4 |grep -v -f - nonref_bubble_coordinate/${genome}.sort.collapse.interectOR.bed0 |cut -f4,8 > nonref_bubble_coordinate/${genome}.vcf_projection
cat nonref_bubble_coordinate/${genome}.vcf_projection | while IFS= read -r line; do
or=$(echo $line |awk '{print $1}')
vb=$(echo $line |awk '{print $2}')
awk -v a=$or '$4==a' anno/${genome}.func.bed > tmp.bed
echo $vb |tr ',' '\n' |grep -f - nonref_bubble_coordinate/${genome}.sort.bed0 |bedtools intersect -b tmp.bed -a - -wb |awk '($3-$2)/($7-$6)==1 {print $8,$4}' 
rm tmp.bed
done > nonref_bubble_coordinate/${genome}.vcf_projection.within
done
```

### Extract non-reference OR's bubbles
```
cat nonref_bubble_coordinate/*.vcf_projection.within |awk '{print $2}' | sort |uniq | while IFS= read -r line; do
echo > subject.fa
count=0
awk -v a=$line '$3==a {print $4","$5}' all.sv.vcf0.vg.correction.vcf0 |tr ',' '\n' | while IFS= read -r line; do
echo -e ">GT$count\n$line" >> subject.fa
count=$(($count+1))
done
minimap2 -cx asm20 all.coding.cdhit98.rename.fasta subject.fa  > aln.paf 
ln=$(awk 'BEGIN{OFS="\t"} {print $1,$3,$4}' aln.paf |sort -k1,1 -k2,2n |bedtools merge -i - -d 30 |awk '$3-$2>800' |wc -l)
if [ $ln -gt 0 ]; then
echo $line
fi
rm subject.fa  aln.paf 
done > nonref_bubble.OR
```

### Genotyping OR for each genotype of all OR bubbles
```
at ref_bubble.OR nonref_bubble.OR | while IFS= read -r line; do
echo > subject.fa
count=0
awk -v a=$line '$3==a {print $4","$5}' all.sv.vcf0.vg.correction.vcf0 |tr ',' '\n' | while IFS= read -r line; do
echo -e ">GT$count\n$line" >> subject.fa
count=$(($count+1))
done
minimap2 -cx asm20 all.coding.cdhit98.rename.fasta subject.fa  > aln.paf 
awk 'BEGIN{OFS="\t"} {print $1,$3,$4,$6}' aln.paf |sort -k1,1 -k2,2n -k4,4 |bedtools merge -i - -c 4 -o collapse -delim "|" |awk '$3-$2>700' > aln.gt
grep '^>' subject.fa |awk -F'>' '{print $2}' | while IFS= read -r gt; do
vb=$(awk -v a=$gt '$1==a' aln.gt)
if [ -n "$vb" ]; then 
awk -v a=$gt '$1==a' aln.gt |cut -f4 |sed ':a;N;$!ba;s/\n/),(/g' |awk -v a=$gt -v b=$line 'BEGIN{OFS="\t"} {print b,a,"("$0")"}'
else
echo -e "$line\t$gt\t*"
fi
done
done > all_bubble.ORgt
```

#### Summarizing OR bubbles' stat
```
for sv in `cut -f1 all_bubble.ORgt |uniq ` ; do
awk -v a=$sv '$1==a' all_bubble.ORgt |cut -f1,3 |sort |uniq |wc -l
done

for sv in `cut -f1 all_bubble.ORgt |uniq ` ; do
tp=$(awk -v a=$sv '$1==a' all_bubble.ORgt |cut -f1,3 |wc -l)
dv=$(awk -v a=$sv '$1==a' all_bubble.ORgt |cut -f3 |sort |uniq|wc -l)
ct=$(awk -v a=$sv '$1==a' all_bubble.ORgt |cut -f3 |sort |uniq |while IFS= read -r line; do
IFS=',' read -r -a elements <<< "$line"
element_count=${#elements[@]}
echo $element_count
done |sort -k1,1 |tail -n 1)
echo -e "$sv\t$tp\t$dv\t$ct"
done > all_bubble.ORgt.stat
```
#### Draw using ggplot2 in R
```
library(ggplot2)
library(tidyr)
df=read.table("all_bubble.ORgt.stat",head=F)
df0=matrix(0,ncol=max(df$V3),nrow=max(df$V2))
for (i in 1:nrow(df)){
a=df[i,2]
b=df[i,3]
df0[a,b]=df0[a,b]+1
}
df0=as.data.frame(df0)
df0$Index <- 1:nrow(df0)
df0_long <- gather(df0, key = "Variable", value = "Value", -Index)
color_scale <- scale_fill_manual(values = colorRampPalette(c("lightblue", "darkblue"))(8))
p <- ggplot(df0_long, aes(x = factor(Index), y = Value, fill = Variable)) +
  geom_bar(stat = "identity") +
  labs(title = "OR bubble genotyping histogram colored with OR discrepency", x = "Number of alternative genotype", y = "Counts", fill = "OR discrepency") +
  theme_minimal() +
  color_scale
ggsave("OR_bubble_hist1.png", plot = p, width = 10, height = 10)

df0=matrix(0,ncol=max(df$V4),nrow=max(df$V2))
for (i in 1:nrow(df)){
a=df[i,2]
b=df[i,4]
df0[a,b]=df0[a,b]+1
}
df0=as.data.frame(df0)
df0$Index <- 1:nrow(df0)
df0_long <- gather(df0, key = "Variable", value = "Value", -Index)
color_scale <- scale_fill_manual(values = colorRampPalette(c("lightblue", "darkblue"))(9))
p <- ggplot(df0_long, aes(x = factor(Index), y = Value, fill = Variable)) +
  geom_bar(stat = "identity") +
  labs(title = "OR bubble genotyping histogram colored with number of ORs", x = "Number of alternative genotype", y = "Counts", fill = "Number of ORs") +
  theme_minimal() +
  color_scale
ggsave("OR_bubble_hist2.png", plot = p, width = 10, height = 10)
```

#### Evaluation of the accuracy of SV genotype by giraffle + call 
GCA_009493645.1 and GCA_009493655.1 assemblies are from the same individual, its NGS reads are aligned and called with sv, theorically NGS's genotype should be same with GCA_009493645.1 + GCA_009493655.1 graph-based genotype

Random selected genotype is used as baseline of accuracy assessment.
```
cut -f1 all_bubble.ORgt |uniq | while IFS= read -r vcfid; do
obsv=$(awk -v a=$vcfid '$3==a' chunk_*.vcf |cut -f10 |awk -F':' '{print $1}')
ob1=$(echo $obsv |awk -F"/" '{print "GT"$1}')
ob2=$(echo $obsv |awk -F"/" '{print "GT"$2}')
rlsv=$(awk -v a=$vcfid '$3==a' all.sv.vcf0.vg.correction.vcf0 |cut -f14,39 | tr '\t' '/')
rl1=$(echo $rlsv |awk -F"/" '{print "GT"$1}')
rl2=$(echo $rlsv |awk -F"/" '{print "GT"$2}')
if [[ "$rl1" != "GT." && "$rl2" != "GT." ]]; then
sh ./CRFD.sh $rl1 $rl2 $ob1 $ob2 $vcfid
else
rl1=$(echo $rlsv |tr '/' '\n'|grep -v "\." |awk '{print "GT"$1}' |head -1)
sh ./CRFD.sh $rl1 $rl1 $ob1 $ob1 $vcfid > a1
sh ./CRFD.sh $rl1 $rl1 $ob2 $ob2 $vcfid > a2
cat a1 a2 | sort -k2,2n |tail -1
fi
done > real_ngs_RC_FD.txt

cut -f1 all_bubble.ORgt |uniq | while IFS= read -r vcfid; do
#obsv=$(awk -v a=$vcfid '$3==a' chunk_*.vcf |cut -f10 |awk -F':' '{print $1}')
ob1=$(awk -v a=$vcfid '$1==a' all_bubble.ORgt |cut -f2 |shuf |head -1)
ob2=$(awk -v a=$vcfid '$1==a' all_bubble.ORgt |cut -f2 |shuf |head -1)
rlsv=$(awk -v a=$vcfid '$3==a' all.sv.vcf0.vg.correction.vcf0 |cut -f14,39 | tr '\t' '/')
rl1=$(echo $rlsv |awk -F"/" '{print "GT"$1}')
rl2=$(echo $rlsv |awk -F"/" '{print "GT"$2}')
if [[ "$rl1" != "GT." && "$rl2" != "GT." ]]; then
sh ./CRFD.sh $rl1 $rl2 $ob1 $ob2 $vcfid
else
rl1=$(echo $rlsv |tr '/' '\n'|grep -v "\." |awk '{print "GT"$1}' |head -1)
sh ./CRFD.sh $rl1 $rl1 $ob1 $ob1 $vcfid > a1
sh ./CRFD.sh $rl1 $rl1 $ob2 $ob2 $vcfid > a2
cat a1 a2 | sort -k2,2n |tail -1
fi
done > random1_ngs_RC_FD.txt
```
