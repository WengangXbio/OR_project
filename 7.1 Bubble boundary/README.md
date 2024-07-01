## Extract bubble's boundary sequence
### Add OR's ID and bubble's length in "all_bubble.ORgt"
```
cut -f1 all_bubble.ORgt |uniq |while IFS= read -r vcfid; do
awk -v a=$vcfid '$1==a'  all_bubble.ORgt |cut -f3 |sort |uniq |while IFS= read -r gt; do
id=$(tr -dc A-Za-z0-9 </dev/urandom | head -c 5)
awk -v a=$vcfid -v b=$gt -v c=$id 'BEGIN{OFS="\t"} $1==a && $3==b {print $0,c}' all_bubble.ORgt
done
done > all_bubble.ORgt.withID

cut -f1 all_bubble.ORgt.withID |uniq |while IFS= read -r vcfid; do
awk -v a=$vcfid '$1==a' all_bubble.ORgt.withID > all_bubble.ORgt.withID.temp
awk -v a=$vcfid '$3==a {print $4}' all.sv.vcf0.vg.correction.vcf0 |awk '{print "GT0",$0}' > gtseq
awk -v a=$vcfid '$3==a {print $5}' all.sv.vcf0.vg.correction.vcf0 |tr ',' '\n' |awk '{print "GT"NR,$0}' >> gtseq
cat all_bubble.ORgt.withID.temp |while IFS= read -r ksk; do
got=$(echo $ksk |awk '{print $2}') 
len=$(awk -v a=$got '$1==a {print $2}' gtseq |wc -m)
echo -e "$ksk\t$len"
done
rm all_bubble.ORgt.withID.temp
done > all_bubble.ORgt.withID.withLen
```

### Extract OR bubble's syntenic regions across 31 cattle assemblies
The strategy is：

1. Looking for syntenic regions of bubbles for each 31 cattle genomes based on MC graph, if there are satisified projections, take the coordinates;
2. If they are bad projections, take minigraph projection as coordinates;
3. Check bubbles's length, which should equal to projection range as above, if the projection range and bubbles length have a difference over 500 bp and difference/length(bubble) >0.1, discard this projection
```
echo "" > Problematic.txt
cut -f1 all_bubble.ORgt.withID.withLen |uniq | while IFS= read -r vcfid; do
	#vcfid=NC_037342.1_78376718
	echo $vcfid
chr=$(awk -v a=$vcfid '$3==a {print $1}' all.sv.vcf0.vg.correction.vcf0)
pos=$(awk -v a=$vcfid '$3==a' all.sv.vcf0.vg.correction.vcf0 |awk '{print $2}')
length=$(awk -v a=$vcfid '$3==a' all.sv.vcf0.vg.correction.vcf0 |awk '{print $4}' |wc -m |awk '{print $1}')
awk -v a=$vcfid '$3==a' all.sv.vcf0.vg.correction.vcf0 |cut -f10- |tr '\t' '\n' |paste genome_prefix - > gtass.tmp
gtind=$(grep -v GCF_002263795 gtass.tmp | awk '$2!="." {print $1}' )
for ind in `echo $gtind` ; do
upid=$(echo "$vcfid"_up"")
downid=$(echo "$vcfid"_down"")
upl=$(awk -v a=$upid '$2==a' ${chr}.untangle.tsv |grep $ind |wc -l)
dwl=$(awk -v a=$downid '$2==a' ${chr}.untangle.tsv |grep $ind |wc -l)
upchr=$(awk -v a=$upid '$2==a' ${chr}.untangle.tsv |grep $ind |awk -F'#' '{print $3}')
dwchr=$(awk -v a=$downid '$2==a' ${chr}.untangle.tsv |grep $ind |awk -F'#' '{print $3}')
upstd=$(awk -v a=$upid '$2==a' ${chr}.untangle.tsv |grep $ind |awk '{print $5}')
dwstd=$(awk -v a=$downid '$2==a' ${chr}.untangle.tsv |grep $ind |awk '{print $5}')
if [[ "$upl" != 1 || "$dwl" != 1 || "$upchr" != "$dwchr" || "$upstd" != "$dwstd" ]] ; then
adj=$(awk -v a=$pos '$2==a' ../1.minigraph/${chr}/${ind}.bed |awk -F'[_:]' '{print $(NF-2)}')
beg=$(awk -v a=$pos '$2==a' ../1.minigraph/${chr}/${ind}.bed |awk -F'[_:]' '{print $(NF-1)}')
end=$(awk -v a=$pos '$2==a' ../1.minigraph/${chr}/${ind}.bed |awk -F'[_:]' '{print $NF}')
beg=$(($beg + $adj - 2000))
end=$(($end + $adj + 2000))
nonchr=$(awk -v a=$pos '$2==a' ../1.minigraph/${chr}/${ind}.bed |awk -F':' '{print $4}' |sed 's/_[^_]*$//')
strand=$(awk -v a=$pos '$2==a' ../1.minigraph/${chr}/${ind}.bed |awk -F':' '{print $3}')
orid=$(awk -v a=$vcfid -v b=$gvd '$1==a && $2==b {print $4}' all_bubble.ORgt.withID.withLen)
if [ "$strand" == "-" ]; then
strand=0
elif [ "$strand" == "+" ]; then
strand=1
fi
gvd=$(grep $ind gtass.tmp |awk '{print "GT"$2}')
aa=$(awk -v a=$gvd -v b=$vcfid '$1==b && $2==a {print $5}' all_bubble.ORgt.withID.withLen)
echo -e "$nonchr\t$beg\t$end\t$ind\t$strand\t$orid\t$gvd\t$aa\tminigraph"
else
nonchr=$(awk -v a=$upid '$2==a' ${chr}.untangle.tsv |grep $ind |awk -F'#' '{print $3}')
adj=$(awk -v a=$upid '$2==a' ${chr}.untangle.tsv |grep $ind |awk -F'#' '{print $4}'|awk -F'[:-]' '{print $2}')
gvd=$(grep $ind gtass.tmp |awk '{print "GT"$2}')
aa=$(awk -v a=$gvd -v b=$vcfid '$1==b && $2==a {print $5}' all_bubble.ORgt.withID.withLen)
orid=$(awk -v a=$vcfid -v b=$gvd '$1==a && $2==b {print $4}' all_bubble.ORgt.withID.withLen)
beg=$(awk -v a=$upid -v b=$downid '$2==a || $2==b' ${chr}.untangle.tsv |grep $ind |cut -f3,4 |tr '\t' '\n' |sort -k1,1n |head -1)
end=$(awk -v a=$upid -v b=$downid '$2==a || $2==b' ${chr}.untangle.tsv |grep $ind |cut -f3,4 |tr '\t' '\n' |sort -k1,1n |tail -n1)
beg=$(($beg + $adj))
end=$(($end + $adj))
echo -e "$nonchr\t$beg\t$end\t$ind\t$upstd\t$orid\t$gvd\t$aa\tMC"
fi
done  > projection.bed
up=$(($pos-2000))
down=$(($pos+$length+2000))
orid=$(awk -v a=$vcfid '$1==a && $2=="GT0" {print $4}' all_bubble.ORgt.withID.withLen)
aa=$(awk -v b=$vcfid '$1==b && $2=="GT0" {print $5}' all_bubble.ORgt.withID.withLen)
echo -e "$chr\t$up\t$down\tGCF_002263795\t1\t$orid\tGT0\t$aa\treference" |cat - projection.bed |sort -k6,6 > projection.sort.bed
rm projection.bed 

awk 'BEGIN{OFS="\t"} {print ($8 > ($3-$2-4000)) ? ($8 - ($3-$2-4000) ) : (($3-$2-4000) - $8), $0}' projection.sort.bed |awk '$1>500 && $1/$9 >0.1' |cut -f2- > problem
cat problem >> Problematic.txt
grep -v -f problem projection.sort.bed > projection.sort.filter.bed
for ind in `awk '$9=="MC" {print $4}' problem` ; do
adj=$(awk -v a=$pos '$2==a' ../1.minigraph/${chr}/${ind}.bed |awk -F'[_:]' '{print $(NF-2)}')
beg=$(awk -v a=$pos '$2==a' ../1.minigraph/${chr}/${ind}.bed |awk -F'[_:]' '{print $(NF-1)}')
end=$(awk -v a=$pos '$2==a' ../1.minigraph/${chr}/${ind}.bed |awk -F'[_:]' '{print $NF}')
beg=$(($beg + $adj - 2000))
end=$(($end + $adj + 2000))
nonchr=$(awk -v a=$pos '$2==a' ../1.minigraph/${chr}/${ind}.bed |awk -F':' '{print $4}' |sed 's/_[^_]*$//')
strand=$(awk -v a=$pos '$2==a' ../1.minigraph/${chr}/${ind}.bed |awk -F':' '{print $3}')
orid=$(awk -v a=$vcfid -v b=$gvd '$1==a && $2==b {print $4}' all_bubble.ORgt.withID.withLen)
if [ "$strand" == "-" ]; then
strand=0
elif [ "$strand" == "+" ]; then
strand=1
fi
gvd=$(grep $ind gtass.tmp |awk '{print "GT"$2}')
aa=$(awk -v a=$gvd -v b=$vcfid '$1==b && $2==a {print $5}' all_bubble.ORgt.withID.withLen)
echo -e "$nonchr\t$beg\t$end\t$ind\t$strand\t$orid\t$gvd\t$aa\tminigraph"
done |awk 'BEGIN{OFS="\t"} {print ($8 > ($3-$2-4000)) ? ($8 - ($3-$2-4000) ) : (($3-$2-4000) - $8), $0}' |awk '$1<500 || $1/$9 <0.1'  |cut  -f2- >> projection.sort.filter.bed
rm problem gtass.tmp
sort -k6,6 projection.sort.filter.bed > projection.sort.filter.bed0

echo "" > boundary_2000_sequence/${vcfid}.sequence.fasta
poly=$(cut -f6 projection.sort.filter.bed0 |sort |uniq |wc -l)
if [ $poly -lt 2 ]; then
echo "Problematic: less_polymorphism: $vcfid" >> Problematic.txt
else
cat projection.sort.filter.bed0 |while IFS= read -r line; do
genome=$(echo $line |awk '{print $4}')
ORID=$(echo $line |awk '{print $6}')
md=$(echo $line |awk '{print $9}')
echo $line| tr ' ' '\t' |bedtools getfasta -fi /public/home/qtyao/WZ_data_transfer/minigraph/0.assembly/${genome}.fasta  -bed - -nameOnly > fasta.tmp
genomerp=$(echo "${genome}_${ORID}_${md}")
sed -i "s/${genome}/${genomerp}/g"  fasta.tmp
strand=$(echo $line |awk '{print $5}')
if [ "$strand" == 0 ]; then
	~/WZ_data_transfer/tools/seqtk-master/seqtk seq -r fasta.tmp 
else
	cat fasta.tmp
fi
rm fasta.tmp
done  >> boundary_2000_sequence/${vcfid}.sequence.fasta
rm  projection.sort.bed projection.sort.filter.bed0
fi
done
```
