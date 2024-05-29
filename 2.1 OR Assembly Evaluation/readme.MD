## Evaluation of OR regions assembly 
### 1. Total number of OR gene in assemblies
### 2. Number of Contigs/scaffold where include OR genes
### 3. Number of OR genes in each chromosome
### 4. Number of Contigs/scaffold composite of each chromosome

#### OR counts and distribution in each chromosome
```
for genome in `ls *.fna`  ; do
numb=$(cut -f1 ${genome}.g2or_itera2_func.bed | uniq |wc -l)
numa=$(cat ${genome}.g2or_itera2_func.bed |wc -l)
echo -e "$numa\t$numb\t$genome"
done
```

#### Count OR by chromosome 
1. processing OR bed annotation by adjusting coordinates
```
cd /share/home/yzwl_zhangwg/OR_project/PGGB/pan_by_chr_new_new
for nchr in `seq 1 30`  ; do
CHR=chr${nchr}
ORdir="/share/home/yzwl_zhangwg/OR_project/genome2OR"
cut -f1 /share/home/yzwl_zhangwg/OR_project/PGGB/pan_by_chr_new_new/${CHR}/${CHR}.pggb.fa.fai > ${CHR}.seq.pathsL
while read line; do
genome=$(echo ${line} | awk -F'[#]' '{print $1}')
anno=$(ls ${ORdir}/*.g2or_itera2_func.bed|grep -v dead |grep $genome)
seq=$(echo ${line} | awk -F'[#]' '{print $3}' |awk -F'[_]' '{print $1}' )
if [ "$seq" == NC ]; then
seq=$(echo ${line} | awk -F'[#]' '{print $3}' |awk -F'[_]' '{print $1"_"$2}' )
start=$(echo ${line} | awk -F'[#]' '{print $3}' |awk -F'[_]' '{print $3}')
end=$(echo ${line} | awk -F'[#]' '{print $3}' |awk -F'[_]' '{print $4}')
else
seq=$(echo ${line} | awk -F'[#]' '{print $3}' |awk -F'[_]' '{print $1}' )
start=$(echo ${line} | awk -F'[#]' '{print $3}' |awk -F'[_]' '{print $2}')
end=$(echo ${line} | awk -F'[#]' '{print $3}' |awk -F'[_]' '{print $3}')
fi
if [ -z "$end" ]; then
awk -v a=$seq -v b=$start -v c=$end -v d=$genome 'BEGIN{OFS="\t"}  $1 == a {print d"#1#"$1,$2,$3,d"_"a"_"$2}' $anno
else
awk -v a=$seq -v b=$start -v c=$end -v d=$genome 'BEGIN{OFS="\t"}  $1 == a && $2 >= b && $3 <= c {print d"#1#"$1"_"b"_"c,$2-b,$3-b,d"_"a"_"$2}' $anno
fi
done < ${CHR}.seq.pathsL > ${CHR}.seq.pathsL.bed
done
```
2. Summary OR counts 
```
fnalist=$(ls *.fasta *.fa *.fna |awk -F'.' '{print $1}')
echo $fnalist |sed "s/ /\n/g" > evaluation_OR_count.txt
for nchr in `seq 1 30`  ; do
for fna in `echo $fnalist` ; do
bed=/share/home/yzwl_zhangwg/OR_project/PGGB/pan_by_chr_new_new/chr${nchr}.seq.pathsL.bed
num=$(grep ${fna} ${bed} |wc -l)
echo $num
done | paste evaluation_OR_count.txt - > temp
mv temp evaluation_OR_count.txt
done 
```

#### Number of Contigs/scaffold composite of each chromosome
```
fnalist=$(ls *.fasta *.fa *.fna |awk -F'.' '{print $1}')
echo $fnalist |sed "s/ /\n/g" > evaluation_chr_scaffold.txt
for i in `seq 1 31`  ; do
chrn=chr$i
for fna in `echo $fnalist`  ; do
ls /share/home/yzwl_zhangwg/OR_project/PGGB/pan_by_chr_new_new/$chrn |awk -F'.' '{print $1}' |uniq -c |grep $fna |awk '{print $1}'
done |paste evaluation_chr_scaffold.txt - > evaluation_chr_scaffold
mv evaluation_chr_scaffold evaluation_chr_scaffold.txt
done
```
