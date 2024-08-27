## PGGB dissection
### Try to solve OR copy number and so forth questions

```
###recreate OR bed with new name, and clue file (all_coding_annotation.rename.bed.clue)
for i in `seq 1 30` ; do
awk -F'#' -v a=${i} '{print $0,$1"_r"a}' ./work/chr${i}/coding_annotation.bed | awk 'BEGIN{OFS="\t"} {print $5"_"$2}' > new_work/temp
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


