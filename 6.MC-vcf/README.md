## Dissection MC vcf by crossfire with MC graph

### 1. preperation vcf and OR annotation
```
mkdir 1.refOR_bub 2.refOR_nobub 3.nonrefOR_proj_bub 4.nonrefOR_proj_nobub 5.nonrefOR_noproj
bedtools intersect -wa -a cattle31.raw.vcfbub.vcf -b GCF_002263795.bed -u > cattle31.raw.vcfbub.vcf.OR
awk 'BEGIN{OFS="\t"} {sum=0; for(i=2; i<=32; i++) sum+=$i; print $0, sum}' ../or.PAV.matrix0.collection |awk '$33!=0 {print $1}' |sed '1d' |sed 's/^[^.]*\.//' > all_existing_OR.name
```

### 1.refOR_bub
This is for the OR within vcf bubble
```
ln -s ../cattle31.raw.vcfbub.vcf.OR ./
ln -s ../GCF_002263795.bed ./
cat cattle31.raw.vcfbub.vcf.OR | while IFS= read -r line; do 
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
done > bubbledOR.gt
cut -f5 1.refOR_bub/bubbledOR.gt | sed 's/,/\n/g'  |awk 'NF' > 1.refOR_bub.ORlist
```

### 2.refOR_nobub
This is for the OR out of vcf bubble
```
ln -s ../1.refOR_bub.ORlist ./
ln -s /public/home/qtyao/WZ_data_transfer/cactus/bos/work/GCF_002263795.func.bed ./
ln -s ../GCF_002263795.bed ./
ln -s ../cattle31.raw.vcfbub.vcf.OR ./
grep -v -f 1.refOR_bub.ORlist GCF_002263795.bed  > refOR_nobub.bed
cut -f4 refOR_nobub.bed | while IFS= read -r line; do
grep ${line} GCF_002263795.func.bed > anno
chr=$(cat anno|awk -F "#" '{print $3}')
start=$(cat anno|awk '{print $2}')
end=$(cat anno|awk '{print $3}')
awk -v a=$chr -v b=$start -v c=$end '$1==a && $2>=b && $2<=c' cattle31.raw.vcfbub.vcf.OR |cut -f10- > genotyping
for ((j=1; j<=30; j++)); do
dot=$(awk -v col=$j '{print $col}' genotyping | grep '\.' |wc -l)
nodot=$(awk -v col=$j '{print $col}' genotyping | grep -v '\.' |wc -l)
echo -e "$dot","$nodot"
done | tr '\n' '\t' | awk -v a=${line} 'BEGIN{OFS="\t"} {print a,$0}'
done > refOR_nobub.gt
cut -f1 2.refOR_nobub/refOR_nobub.gt > 2.refOR_nobub.ORlist
```

### 3.nonrefOR_proj_bub
This is for the nonreference OR within vcf bubble
```
awk 'BEGIN{OFS="\t"} {sum=0; for(i=2; i<=32; i++) sum+=$i; print $0, sum}' or.PAV.matrix0.collection |awk '$33!=0 {print $1}' |sed '1d' |sed 's/^[^.]*\.//' > all_varied.ORlist
ln -s ../all_varied.ORlist ./
ln -s ../../../or.PAV.matrix0.collection ./
ln -s ../../../seq.txt ./
ln -s ../cattle31.raw.vcfbub.vcf ./
grep -v "GCF_002263795" all_varied.ORlist > nonref.ORlist

cat nonref.ORlist | while IFS= read -r line; do 
OR=${line}
chr=$(grep ${OR} or.PAV.matrix0.collection |awk -F'.' '{print $1}')
cotig=$(awk -v a=${chr} '$1==a' seq.txt |cut -f2)
grep ${OR} ../../../chr_${cotig}/injref.untangle.tsv.collapse |grep GCF_002263795
done > nonrefOR_proj.injref.untangle.tsv.collapse
awk -F['#:\t'] 'BEGIN{OFS="\t"} {print $3,$6,$7,$5}' nonrefOR_proj.injref.untangle.tsv.collapse |sort -k1,1 -k2,2n > nonrefOR_proj.injref.untangle.tsv.collapse.bed

bedtools intersect -wa -a cattle31.raw.vcfbub.vcf -b nonrefOR_proj.injref.untangle.tsv.collapse.bed -u > cattle31.raw.vcfbub.vcf.nonrefOR_proj
cat cattle31.raw.vcfbub.vcf.nonrefOR_proj | while IFS= read -r line; do 
chr=$(echo $line |awk '{print $1}')
str=$(echo $line |awk '{print $2}')
len=$(echo $line |awk '{print $4}' |awk '{print length($0)}')
end=$(($str+$len)); alt=$(echo $line |awk '{print $5}')
alt_len=$(echo $alt | awk -F, '{for(i=1;i<=NF;i++) printf "%d%s", length($i), (i<NF ? "," : "\n")}')
cor_or=$(echo -e "$chr\t$str\t$end" | bedtools intersect -b nonrefOR_proj.injref.untangle.tsv.collapse.bed -a - -wb |awk '($3-$2)/($6-$5) ==1 {print $7}' |tr '\n' ',')
gt=$(echo $line |cut -d' ' -f10-)
if [ -n "$cor_or" ]; then 
echo -e "$chr\t$str\t$len\t$alt_len\t$cor_or\t$gt"
fi
done > nonrefOR_proj.bubbledOR.gt
cut -f5 3.nonrefOR_proj_bub/nonrefOR_proj.bubbledOR.gt | sed 's/,/\n/g'  |awk 'NF' > 3.nonrefOR_proj_bub.ORlist
```

### 4.nonrefOR_proj_nobub
This is for the nonreference OR outof vcf bubble
```
ln -s ../3.nonrefOR_proj_bub.ORlist ./
ln -s ../3.nonrefOR_proj_bub/nonrefOR_proj.injref.untangle.tsv.collapse.bed ./
ln -s ../3.nonrefOR_proj_bub/cattle31.raw.vcfbub.vcf.nonrefOR_proj ./
grep -v -f 3.nonrefOR_proj_bub.ORlist nonrefOR_proj.injref.untangle.tsv.collapse.bed  > nonrefOR_proj_nobub.bed
cat nonrefOR_proj_nobub.bed | while IFS= read -r line; do
chr=$(echo ${line}|awk '{print $1}')
start=$(echo ${line}|awk '{print $2}')
end=$(echo ${line}|awk '{print $3}')
or=$(echo ${line}|awk '{print $4}')
awk -v a=$chr -v b=$start -v c=$end '$1==a && $2>=b && $2<=c' cattle31.raw.vcfbub.vcf.nonrefOR_proj |cut -f10- > genotyping
for ((j=1; j<=30; j++)); do
dot=$(awk -v col=$j '{print $col}' genotyping | grep '\.' |wc -l)
nodot=$(awk -v col=$j '{print $col}' genotyping | grep -v '\.' |wc -l)
echo -e "$dot","$nodot"
done | tr '\n' '\t' | awk -v a=${or} 'BEGIN{OFS="\t"} {print a,$0}'
done > nonrefOR_proj_nobub.gt
cut -f1 4.nonrefOR_proj_nobub/nonrefOR_proj_nobub.gt > 4.nonrefOR_proj_nobub.ORlist
```

### 5.nonrefOR_noproj
```
ln -s ../3.nonrefOR_proj_bub/nonref.ORlist ./
ln -s ../3.nonrefOR_proj_bub.ORlist ./
ln -s ../4.nonrefOR_proj_nobub.ORlist ./
ln -s ../../../all.func.bed ./
ln -s ../../../../genomelist ./
ln -s ../../../or.PAV.matrix0.collection ./
ln -s ../../../seq.txt ./
ln -s ../cattle31.raw.vcfbub.vcf ./
cat 3.nonrefOR_proj_bub.ORlist 4.nonrefOR_proj_nobub.ORlist |grep -v -f - nonref.ORlist > nonrefOR_noproj.ORlist

cat nonrefOR_noproj.ORlist | while IFS= read -r OR; do 
genome=$(cat all.func.bed |grep $OR |awk -F'#' '{print $1}')
chr=$(cat all.func.bed |grep $OR |awk -F'#' '{print $3}')
start=$(($(cat all.func.bed |grep $OR |awk '{print $2}') - 50000))
end=$(($(cat all.func.bed |grep $OR |awk '{print $3}')+ 50000))
if [ "$start" -lt 0 ]; then 
start=0
fi
fasta=$(grep $genome genomelist |awk '{print $2}')
echo -e "$chr\t$start\t$end\t$OR" |bedtools getfasta -name -fi ${fasta} -bed - 
done > nonrefOR_noproj.ORlist.fasta
fasta=/public/home/qtyao/WZ_data_transfer/cactus/bos/assembly/GCF_002263795.3_ARS-UCD2.0_genomic.fna
minimap2 -cx asm5 ${fasta} nonrefOR_noproj.ORlist.fasta  -t 4 >  nonrefOR_noproj.ORlist.fasta.paf

for OR in `cat nonrefOR_noproj.ORlist ` ; do
chr=$(grep ${OR} or.PAV.matrix0.collection |awk -F'.' '{print $1}')
cotig=$(awk -v a=${chr} '$1==a' seq.txt |cut -f2)
coor=$(awk -v a=$OR -v b=$cotig '$1==a && $6==b' nonrefOR_noproj.ORlist.fasta.paf |cut -f8,9 |tr '\t' '\n' |sort -k1,1n)
while [ -n "$coor" ]; do 
start=$(echo $coor |tr ' ' '\n' |head -1)
end=$(echo $coor |tr ' ' '\n' |awk -v a=$start '$1<a+220000' |tail -n1)
echo -e "$cotig\t$start\t$end\t$OR" 
coor=$(echo $coor |tr ' ' '\n' |awk -v a=$start '$1>a+220000')
done
done | awk 'BEGIN{OFS="\t"} {print $1,$2-10000,$3+10000,$4}' > nonrefOR_noproj.ORlist.projected_ref.bed

cat nonrefOR_noproj.ORlist.projected_ref.bed | while IFS= read -r line; do
chr=$(echo $line |awk '{print $1}')
start=$(echo $line |awk '{print $2}')
end=$(echo $line |awk '{print $3}')
OR=$(echo $line |awk '{print $4}')
bsub -J insertion_bubble_finder -n 1 -R "span[hosts=1] rusage[mem=2GB]" -o %J.out -e %J.err -q normal \
"
sh ./insertion_bubble_finder.sh ${chr} $start $end $OR
"
done 
cat *_nonrefOR_noproj.insertion_bubble > nonrefOR_noproj.insertion_bubble

cut -f1,2 nonrefOR_noproj.insertion_bubble |sort |uniq | while IFS= read -r line; do 
chr=$(echo $line |awk '{print $1}')
pos=$(echo $line |awk '{print $2}')
vcf=$(awk -v a=${chr} -v b=${pos} '$1==a && $2==b' cattle31.raw.vcfbub.vcf)
len=$(echo $vcf |awk '{print $4}' |awk '{print length($0)}')
alt_len=$(echo $vcf |awk '{print $5}'  | awk -F, '{for(i=1;i<=NF;i++) printf "%d%s", length($i), (i<NF ? "," : "\n")}')
gt=$(echo $vcf |cut -d' ' -f10-)
or=$(awk -v a=${chr} -v b=${pos} '$1==a && $2==b' nonrefOR_noproj.insertion_bubble |cut -f3 |tr '\n' ',')
echo -e "$chr\t$pos\t$len\t$alt_len\t$or\t$gt"
done > nonrefOR_noproj.bubbledOR.gt

cut -f3 5.nonrefOR_noproj/nonrefOR_noproj.insertion_bubble > 5.nonrefOR_noproj.ORlist
```
### 6.all_collection
#### 6.1 correction graph projection matrix by rescan inject.untangle.tsv
```
cat bubbledOR.gt nonrefOR_noproj.bubbledOR.gt nonrefOR_proj.bubbledOR.gt > all_bubbledOR.gt
cat nonrefOR_proj_nobub.gt refOR_nobub.gt > all_nobubbledOR.gt

sed '1d' or.PAV.matrix0.collection | while IFS= read -r line; do 
col1=$(echo $line |awk '{print $1}')
chr=$(echo $line |awk -F'.' '{print $1}')
or=$(echo $line |awk '{print $1}' |sed 's/^[^.]*\.//')
ass=$(awk -v a=${chr} '$1==a {print $2}' seq.txt)
graphgt=$(echo $line|cut -d' ' -f2-)
echo $graphgt |tr ' ' '\n'|paste or.PAV.matrix0.collection.index - > tegraphgt
head -n31 tegraphgt | while IFS= read -r tgg ; do 
genome=$(echo $tgg |awk '{print $1}' )
gt=$(echo $tgg |awk '{print $2}')
if  [ $gt == 0 ]; then 
len=$(awk -v a=${or} '$2==a' ../../../chr_${ass}/injref.untangle.tsv |grep ^${genome}# |uniq| awk '{print $4-$3}'| awk '{ sum += $1} END {print sum}' |awk '$1>=100')
ln=$(awk -v a=${or} '$2==a' ../../../chr_${ass}/injref.untangle.tsv |grep ^${genome}# |uniq| awk '{print $4-$3}'| awk '{ sum += $1} END {print sum}' |awk '$1>=100' |wc -l)
echo $ln
if  [ $gt != $ln ]; then 
echo -e "$or\t$genome\t$gt\t$ln\t$len" >> or.PAV.matrix0.collection.new.log
fi
else
echo $gt
fi
done | tr '\n' '\t' |awk -v a=$col1 'BEGIN{OFS="\t"} {print a,$0}'
done > or.PAV.matrix0.collection.new
head -n1 or.PAV.matrix0.collection |cat - or.PAV.matrix0.collection.new > or.PAV.matrix0.collection.new0
```
#### 6.2 Detection of conflict between graph projection and vcf genotyping
(1) Bubble vcf genotyping

(2) No bubble vcf genotyping

(3) Same bubble vcf genotype with each other
```
cut -f1 all_nobubbledOR.gt |grep -f - or.PAV.matrix0.collection | while IFS= read -r line; do 
chr=$(echo $line |awk -F'.' '{print $1}')
or=$(echo $line |awk '{print $1}' |sed 's/^[^.]*\.//')
vcfgt=$(grep $or all_nobubbledOR.gt|cut -f2-)
graphgt=$(echo $line|cut -d' ' -f2-)
graphgt0=$(echo $graphgt |tr ' ' '\n' |paste or.PAV.matrix0.collection.index - |grep -v GCF_002263795 |awk '$2==0 {print $1}')
echo $vcfgt | tr ' ' '\n' | paste gt.index - > vcftype
if [ -n "$graphgt0" ]; then 
echo $graphgt0 | tr ' ' '\n' |grep -f -  vcftype | awk -v b=${or} -v a=${chr} -F'[\t,]' '$2<$3 {print a,b,$0}' | while IFS= read -r line; do 
chr=$(echo $line |awk '{print $1}')
ass=$(awk -v a=${chr} '$1==a {print $2}' seq.txt)
or=$(echo $line |awk '{print $2}')
genome=$(echo $line |awk '{print $3}')
ln=$(grep ${or} ../../../chr_${ass}/injref.untangle.tsv |grep ${genome} |wc -l)
if [ $ln == 0 ]; then
echo $line
fi
done
fi
done > all_nobubbledOR.conflict


cut -f5 all_bubbledOR.gt |tr ',' '\n' |grep -v "^$" |grep -f - or.PAV.matrix0.collection.new0 | while IFS= read -r line; do 
chr=$(echo $line |awk -F'.' '{print $1}')
or=$(echo $line |awk '{print $1}' |sed 's/^[^.]*\.//')
vcfgt=$(grep $or all_bubbledOR.gt|cut -f6-)
echo $vcfgt | tr ' ' '\n' |awk '{print $1}' | paste gt.index - > vcftype
vcfgtlt=$(grep $or all_bubbledOR.gt|cut -f4 |tr ',' '\n' |awk '{print NR,$0}' |awk '$2<100 {print $1}' )
graphgt=$(echo $line|cut -d' ' -f2-)
for gtt in `echo -e "$vcfgtlt\n".""` ; do
awk -v a=${gtt} '$2==a {print $1}' vcftype
done > vcfgtsid
echo $graphgt |tr ' ' '\n' |paste or.PAV.matrix0.collection.new0.index - |grep -f vcfgtsid - |awk -v a=$or '$2!=0 {print a,$0}' > confliction
awk '{print $2}' confliction |grep -f - vcftype |awk '{print $2}' |paste confliction -
done > all_bubbledOR.conflict


cut -f5 all_bubbledOR.gt |tr ',' '\n' |grep -v "^$" |grep -f - or.PAV.matrix0.collection.new0 | while IFS= read -r line; do 
chr=$(echo $line |awk -F'.' '{print $1}')
or=$(echo $line |awk '{print $1}' |sed 's/^[^.]*\.//')
vcfgt=$(grep $or all_bubbledOR.gt|cut -f6-)
echo $vcfgt | tr ' ' '\n' |awk '{print $1}' | paste gt.index - > vcftype
congt=$(awk '{print $2}' vcftype  |sort |uniq -c |grep -v "\." |awk '$1>1 {print $2}')
graphgt=$(echo $line|cut -d' ' -f2-)
for lds in `echo $congt` ; do
awk -v a=${lds} '$2==a {print $1}' vcftype >ldsid
con=$(echo $graphgt |tr ' ' '\n' |paste or.PAV.matrix0.collection.new0.index - | grep -f ldsid - |awk '{print $2}' |sort |uniq |wc -l)
if [ $con != 1 ]; then
echo -e "$or\t$lds"
fi
done
done > all_bubbledOR.concensus.conflict
```

#### 6.3 Create vcf-based matrix
```
sed '1d' or.PAV.matrix0.collection.new0 | while IFS= read -r line; do 
wd=$(echo $line |awk '{print $1}')
or=$(echo $line |awk '{print $1}' |sed 's/^[^.]*\.//')
tt=$(grep $or all_bubbledOR.gt)
if [ -n "$tt" ]; then
echo $line|cut -d' ' -f2-  |tr ' ' '\n' |paste or.PAV.matrix0.collection.new0.index - > graph.gt
vcfgt=$(grep $or all_bubbledOR.gt|cut -f6-)
echo $vcfgt | tr ' ' '\n' |awk '{print $1}' | paste gt.index - |awk '$2=="." {print $1}' > missid
for id in `cat missid ` ; do
awk -v a=${id} 'BEGIN{OFS="\t"} {if ($1==a) {print $1,"NA"}else {print $0} }' graph.gt > tmp
mv tmp graph.gt
done
cut -f2 graph.gt | tr '\n' '\t' |awk -v a=$wd 'BEGIN{OFS="\t"} {print a,$0}'
fi
tt=$(grep $or all_nobubbledOR.gt)
if [ -n "$tt" ]; then
echo $line|cut -d' ' -f2-  |tr ' ' '\n' |paste or.PAV.matrix0.collection.new0.index - > graph.gt
vcfgt=$(grep $or all_nobubbledOR.gt|cut -f2-)
echo $vcfgt | tr ' ' '\n' |awk '{print $1}' | paste gt.index - |awk -F'[\t,]'  '$2>$3 {print $1}' > missid
for id in `cat missid ` ; do
awk -v a=${id} 'BEGIN{OFS="\t"} {if ($1==a) {print $1,"NA"}else {print $0} }' graph.gt > tmp
mv tmp graph.gt
done
cut -f2 graph.gt | tr '\n' '\t' |awk -v a=$wd 'BEGIN{OFS="\t"} {print a,$0}'
fi
done > or.PAV.matrix0.collection.new0.vcfNA
head -n1 or.PAV.matrix0.collection |cat - or.PAV.matrix0.collection.new0.vcfNA > or.PAV.matrix0.collection.new0.vcfNA0
```

#### 6.4 Create graph-based matrix
```
sed '1d' or.PAV.matrix0.collection.new0 | while IFS= read -r line; do 
wd=$(echo $line |awk '{print $1}')
or=$(echo $line |awk '{print $1}' |sed 's/^[^.]*\.//')
tt=$(grep $or all_bubbledOR.gt)
if [ -n "$tt" ]; then
echo $line|cut -d' ' -f2-  |tr ' ' '\n' |paste or.PAV.matrix0.collection.new0.index - > graph.gt
vcfgt=$(grep $or all_bubbledOR.gt|cut -f6-)
echo $vcfgt | tr ' ' '\n' |awk '{print $1}' | paste gt.index - |awk '$2=="." {print $1}' > missid
for id in `cat missid ` ; do
awk -v a=${id} 'BEGIN{OFS="\t"} {if ($1==a && $2==0) {print $1,"NA"}else {print $0} }' graph.gt > tmp
mv tmp graph.gt
done
cut -f2 graph.gt | tr '\n' '\t' |awk -v a=$wd 'BEGIN{OFS="\t"} {print a,$0}'
fi
tt=$(grep $or all_nobubbledOR.gt)
if [ -n "$tt" ]; then
echo $line|cut -d' ' -f2-  |tr ' ' '\n' |paste or.PAV.matrix0.collection.new0.index - > graph.gt
vcfgt=$(grep $or all_nobubbledOR.gt|cut -f2-)
echo $vcfgt | tr ' ' '\n' |awk '{print $1}' | paste gt.index - |awk -F'[\t,]'  '$2>$3 {print $1}' > missid
for id in `cat missid ` ; do
awk -v a=${id} 'BEGIN{OFS="\t"} {if ($1==a && $2==0) {print $1,"NA"}else {print $0} }' graph.gt > tmp
mv tmp graph.gt
done
cut -f2 graph.gt | tr '\n' '\t' |awk -v a=$wd 'BEGIN{OFS="\t"} {print a,$0}'
fi
done > or.PAV.matrix0.collection.new0.graphNA
head -n1 or.PAV.matrix0.collection |cat - or.PAV.matrix0.collection.new0.graphNA > or.PAV.matrix0.collection.new0.graphNA0
```
#### 6.5 Create vcf-based matrix only for PAV OR
```
sed '1d' or.PAV.matrix0.collection.new0 | while IFS= read -r line; do 
wd=$(echo $line |awk '{print $1}')
or=$(echo $line |awk '{print $1}' |sed 's/^[^.]*\.//')
tt=$(grep $or all_bubbledOR.gt)
if [ -n "$tt" ]; then
echo $line|cut -d' ' -f2-  |tr ' ' '\n' |paste or.PAV.matrix0.collection.new0.index - > graph.gt
vcfgt=$(grep $or all_bubbledOR.gt|cut -f6-)
echo $vcfgt | tr ' ' '\n' |awk '{print $1}' | paste gt.index - |awk '$2=="." {print $1}' > missid
for id in `cat missid ` ; do
awk -v a=${id} 'BEGIN{OFS="\t"} {if ($1==a) {print $1,"NA"}else {print $0} }' graph.gt > tmp
mv tmp graph.gt
done
cut -f2 graph.gt | tr '\n' '\t' |awk -v a=$wd 'BEGIN{OFS="\t"} {print a,$0}'
fi
done > or.PAV.matrix0.collection.new0.vcfNA.PAV
head -n1 or.PAV.matrix0.collection |cat - or.PAV.matrix0.collection.new0.vcfNA.PAV > or.PAV.matrix0.collection.new0.vcfNA.PAV0
```
