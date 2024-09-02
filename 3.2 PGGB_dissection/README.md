## PGGB dissection
#### Try to solve OR copy number and so forth questions

### 1. reannotation
```
###recreate OR bed with new name, and clue file (all_coding_annotation.rename.bed.clue)
cd /share/org/YZWL/yzwl_zhangwg/OR_project/1.pangenome/1.3bos_PGGB/pan_cattle31_new
for i in `seq 1 30` ; do
awk -F'#' -v a=${i} '{print $0,$1"_R"a}' ./work/chr${i}/coding_annotation.bed | awk 'BEGIN{OFS="\t"} {print $5"_"$2}' > new_work/temp
cat new_work/replace_genome.clue |while IFS= read -r line; do
old=$(echo $line |awk '{print $1}')
new=$(echo $line |awk '{print $2}')
sed -i "s/${old}/${new}/g"  new_work/temp
done
awk 'BEGIN{OFS="\t"} {print $4}' ./work/chr${i}/coding_annotation.bed |paste - new_work/temp >> new_work/all_coding_annotation.rename.bed.clue
awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' ./work/chr${i}/coding_annotation.bed |paste - new_work/temp 
done > new_work/all_coding_annotation.rename.bed

###Extarct fasta based on bed
for i in `seq 1 30` ; do
bedtools getfasta -fi chr${i}.pggb.fa -bed new_work/all_coding_annotation.rename.bed -nameOnly
done  > new_work/all_coding_annotation.rename.fasta

###Change lower case into capital
tr 'a-z' 'A-Z' < all_coding_annotation.rename.fasta > all_coding_annotation.rename.FASTA

###Reverse-complement backward sequences
awk '/^ATG/{print prev} {prev=$0}' all_coding_annotation.rename.FASTA |awk -F'>' '{print $2}' > all_coding_annotation.rename.forward
grep "^>" all_coding_annotation.rename.FASTA |awk -F'>' '{print $2}' |grep -v -f all_coding_annotation.rename.forward - > all_coding_annotation.rename.backward
grep -A 1 -wFf all_coding_annotation.rename.forward all_coding_annotation.rename.FASTA |grep -v "-" > all_coding_annotation.rename.forward.FASTA
grep -A 1 -wFf all_coding_annotation.rename.backward all_coding_annotation.rename.FASTA |grep -v "-" > all_coding_annotation.rename.backward.FASTA
/share/org/YZWL/yzwl_zhangwg/tools/seqtk/seqtk seq -r all_coding_annotation.rename.backward.FASTA > all_coding_annotation.rename.backward.rc.FASTA
cat all_coding_annotation.rename.forward.FASTA all_coding_annotation.rename.backward.rc.FASTA > all_coding_annotation.rename.ordered.FASTA
```

### 2. projection to calculate orthology
```
cd /share/org/YZWL/yzwl_zhangwg/OR_project/1.pangenome/1.3bos_PGGB/pan_cattle31_new/new_work/2.projection
for i in `seq 1 30` ; do
index=1
grep -v '^$' chr${i}.group |while IFS= read -r OR; do
printf -v index_padded "%03d" $index
ORid=$(echo R${i}_${index_padded})

####### write coding OR #######
echo $OR |tr ',' '\n'|grep -v '^$' | grep -f - ../all_coding_annotation.rename.bed.clue |awk -v a=${ORid} '$1=a {print $0}' >> coding.anno

####### write validated pseudo OR #######
echo $OR |tr ',' '\n'|grep -v '^$' | \
grep -f - /share/org/YZWL/yzwl_zhangwg/OR_project/1.pangenome/1.3bos_PGGB/pan_cattle31_new/work/chr${i}/injref.untangle.tsv.collapse | \
awk 'BEGIN{OFS="\t"} {print $1,$3,$4}' | bedtools sort -i - |\
bedtools subtract -A -a - -b /share/org/YZWL/yzwl_zhangwg/OR_project/1.pangenome/1.3bos_PGGB/pan_cattle31_new/work/chr${i}/coding_annotation.bed | \
bedtools merge -i - |\
bedtools intersect -wa -b -  -a /share/org/YZWL/yzwl_zhangwg/OR_project/1.pangenome/1.3bos_PGGB/pan_cattle31_new/work/chr${i}/dead_annotation.bed > validated_pseudo_proto.bed
awk -F'#' -v a=$i '{print $0,$1"_R"a}' validated_pseudo_proto.bed | awk 'BEGIN{OFS="\t"} {print $5"_"$2}' > temp 
cat ../replace_genome.clue |while IFS= read -r line; do
old=$(echo $line |awk '{print $1}')
new=$(echo $line |awk '{print $2}')
sed -i "s/${old}/${new}/g"  temp
done
awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' validated_pseudo_proto.bed |paste - temp >> validated_pseudo.bed
awk 'BEGIN{OFS="\t"} {print $4}' validated_pseudo_proto.bed |paste - temp >> validated_pseudo.bed.clue
awk -v a=$ORid 'BEGIN{OFS="\t"} {print a,$1}' temp >> validated_pseudo.bed.anno

####### write unvalidated pseudo OR #######
echo $OR |tr ',' '\n'|grep -v '^$' | \
grep -f - /share/org/YZWL/yzwl_zhangwg/OR_project/1.pangenome/1.3bos_PGGB/pan_cattle31_new/work/chr${i}/injref.untangle.tsv.collapse | \
awk 'BEGIN{OFS="\t"} {print $1,$3,$4}' | bedtools sort -i - |\
bedtools subtract -A -a - -b /share/org/YZWL/yzwl_zhangwg/OR_project/1.pangenome/1.3bos_PGGB/pan_cattle31_new/work/chr${i}/coding_annotation.bed | \
bedtools merge -i - |\
bedtools subtract -A -a -  -b /share/org/YZWL/yzwl_zhangwg/OR_project/1.pangenome/1.3bos_PGGB/pan_cattle31_new/work/chr${i}/dead_annotation.bed  > temp 
awk -F '[#\t]'  -v a=$i -v b=$ORid  'BEGIN{OFS="\t"} {print b,$1"_R"a"_"$4}' temp > temp.anno
cat ../replace_genome.clue |while IFS= read -r line; do
old=$(echo $line |awk '{print $1}')
new=$(echo $line |awk '{print $2}')
sed -i "s/${old}/${new}/g"  temp.anno
done
cut -f2 temp.anno | paste temp - >> unvalidated_pseudo.bed
cat temp.anno >> unvalidated_pseudo.bed.anno

index=$((index + 1))
done
done
```

### 3. Redefine OR orthology group
Some orthology defined by graph is not accurate,  regroup OR by cd-hit cluster, 7 OR are regroup
```
cd /share/org/YZWL/yzwl_zhangwg/OR_project/1.pangenome/1.3bos_PGGB/pan_cattle31_new/new_work/3.classification
```
See "manual_anno.sh"

### 4. OR statistic 
```
#calculate all 31 bos genome (27 cattle + 4 wild cattle)
for OR in `awk '{print $1}'  coding.anno.update |uniq` ; do
all_num=$(cat validated_pseudo.bed.anno unvalidated_pseudo.bed.anno coding.anno.update |awk -v a=$OR '$1==a {print $2}' |wc -l)
coding_num=$(awk -v a=$OR '$1==a {print $2}' coding.anno.update |wc -l)
psudo_num=$(cat validated_pseudo.bed.anno unvalidated_pseudo.bed.anno |awk -v a=$OR '$1==a {print $2}' |wc -l)
coding_genome_num=$(awk -v a=$OR '$1==a {print $2}' coding.anno.update |awk -F'_' '{print $1}' |sort |uniq |wc -l)
psudo_genome_num=$(cat validated_pseudo.bed.anno unvalidated_pseudo.bed.anno | awk  -v a=$OR '$1==a {print $2}' |awk -F'_' '{print $1}' |sort |uniq |wc -l)
all_genome_num=$(cat coding.anno.update validated_pseudo.bed.anno unvalidated_pseudo.bed.anno | awk  -v a=$OR '$1==a {print $2}' |awk -F'_' '{print $1}' |sort |uniq |wc -l)
echo -e "$OR\t$all_num\t$coding_num\t$psudo_num\t$all_genome_num\t$coding_genome_num\t$psudo_genome_num" 
done > coding.anno.update.stat

#calculate 27 cattle genome
echo "A5 A26 A29 A30 A31" |tr ' ' '\n' |grep  -v -f - coding.anno.update > coding.anno.update.cattle
echo "A5 A26 A29 A30 A31" |tr ' ' '\n' |grep  -v -f - validated_pseudo.bed.anno > validated_pseudo.bed.anno.cattle
echo "A5 A26 A29 A30 A31" |tr ' ' '\n' |grep  -v -f - unvalidated_pseudo.bed.anno > unvalidated_pseudo.bed.anno.cattle

for OR in `awk '{print $1}'  coding.anno.update.cattle |uniq` ; do
all_num=$(cat validated_pseudo.bed.anno.cattle unvalidated_pseudo.bed.anno.cattle coding.anno.update.cattle |awk -v a=$OR '$1==a {print $2}' |wc -l)
coding_num=$(awk -v a=$OR '$1==a {print $2}' coding.anno.update.cattle |wc -l)
psudo_num=$(cat validated_pseudo.bed.anno.cattle unvalidated_pseudo.bed.anno.cattle |awk -v a=$OR '$1==a {print $2}' |wc -l)
coding_genome_num=$(awk -v a=$OR '$1==a {print $2}' coding.anno.update.cattle |awk -F'_' '{print $1}' |sort |uniq |wc -l)
psudo_genome_num=$(cat validated_pseudo.bed.anno.cattle unvalidated_pseudo.bed.anno.cattle | awk  -v a=$OR '$1==a {print $2}' |awk -F'_' '{print $1}' |sort |uniq |wc -l)
all_genome_num=$(cat coding.anno.update.cattle validated_pseudo.bed.anno.cattle unvalidated_pseudo.bed.anno.cattle | awk  -v a=$OR '$1==a {print $2}' |awk -F'_' '{print $1}' |sort |uniq |wc -l)
echo -e "$OR\t$all_num\t$coding_num\t$psudo_num\t$all_genome_num\t$coding_genome_num\t$psudo_genome_num" 
done > coding.anno.update.stat.cattle

#single-copy
awk '$2==$5  ' coding.anno.update.stat.cattle |wc -l
#single-copy
```

### 5. Are orthology genes collasped in graph?
```
cd /share/org/YZWL/yzwl_zhangwg/OR_project/1.pangenome/1.3bos_PGGB/pan_cattle31_new/new_work/1.annotation
grep '^GCF_002263795' all_coding_annotation.rename.bed.clue > /share/org/YZWL/yzwl_zhangwg/OR_project/1.pangenome/1.3bos_PGGB/pan_cattle31_new/new_work/5.orthology_collasped_graph/ARS-UCD2.0.graph.clue
cd /share/org/YZWL/yzwl_zhangwg/OR_project/2.outi_evolution/1.annotation/annotation
awk -F'[_\t]' 'BEGIN{OFS="\t"}{print $2,$3,$4,$7}' GCF_002263795.3_ARS-UCD2.0_genomic.fna.func_pro_ORs.fasta.rename.trace |sed 's/@#@/_/g' - > /share/org/YZWL/yzwl_zhangwg/OR_project/1.pangenome/1.3bos_PGGB/pan_cattle31_new/new_work/5.orthology_collasped_graph/ARS-UCD2.0.outi.bed
cd /share/org/YZWL/yzwl_zhangwg/OR_project/1.pangenome/1.3bos_PGGB/pan_cattle31_new/new_work/5.orthology_collasped_graph
sh modify_coor.sh ARS-UCD2.0.outi.bed ARS-UCD2.0.outi.mbed
awk 'BEGIN{OFS="\t"} {print "GCF_002263795_"$1"_"$2,$4}' ARS-UCD2.0.outi.mbed > ARS-UCD2.0.outi.clue
grep "NC_0373" ARS-UCD2.0.outi.clue |sort -k1,1 > ARS-UCD2.0.outi.sort.clue
sort -k1,1 ARS-UCD2.0.graph.clue  > ARS-UCD2.0.graph.sort.clue
paste ARS-UCD2.0.outi.sort.clue ARS-UCD2.0.graph.sort.clue |cut -f1,2,4 > ARS-UCD2.0.graph.outi.link
awk 'BEGIN{OFS="\t"} NR==FNR { map[$2] = $1; next } { print $0, (map[$3] ? map[$3] : "") }' coding.anno.update ARS-UCD2.0.graph.outi.link > ARS-UCD2.0.graph.outi.link.group
```

### 6. OR distance and similarity
```
grep "^>" all_coding_annotation.rename.ordered.FASTA |grep "A27" |awk -F'>' '{print $2}' > A27.OR.name
~/tools/seqtk/seqtk subseq all_coding_annotation.rename.ordered.FASTA A27.OR.name > A27.OR.fasta
python calculate_similarity.py A27.OR.fasta > A27.OR.similarity

awk -F'[_ ]+' '{print $2,$3,$5,$6,$7}' A27.OR.similarity > A27.OR.similarity.col

#R
df=read.table("A27.OR.similarity.col",head=F)
df1=df[which(df$V1==df$V3),]
model =lm(df1$V5~df1$V4-df1$V2)
cor.test(abs(df1$V4-df1$V2),df1$V5)
 plot(abs(df1$V4-df1$V2),df1$V5,xlab="Distance",ylab="Similarity")
abline(model, col="blue",lwd=4)
```
