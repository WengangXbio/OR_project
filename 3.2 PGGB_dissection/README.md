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
