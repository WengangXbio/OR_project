## Extract bubble's boundary sequence
### 1. Extract projection region in graph
```
cut -f1 all_bubble.ORgt |uniq |while IFS= read -r vcfid; do
chr=$(awk -v a=$vcfid '$3==a' all.sv.vcf0.vg.correction.vcf0 |awk '{print $1}')
pos=$(awk -v a=$vcfid '$3==a' all.sv.vcf0.vg.correction.vcf0 |awk '{print $2}')
length=$(awk -v a=$vcfid '$3==a' all.sv.vcf0.vg.correction.vcf0 |awk '{print $4}' |wc -m |awk '{print $1}')
seq=$(grep $chr /public/home/qtyao/WZ_data_transfer/cactus/bos/work/cattle31.${chr}.full.og.L)
up1=$(($pos-2000))
down1=$(($pos+$length))
down2=$(($down1+2000))
echo -e "$seq\t$up1\t$pos\t$vcfid"_up""    >> ${chr}.bed
echo -e "$seq\t$down1\t$down2\t$vcfid"_down""  >> ${chr}.bed
done 


chr="NC_037357 NC_037356 NC_037355 NC_037352 NC_037351 NC_037350 NC_037346 NC_037343 NC_037342 NC_037338 NC_037337 NC_037336 NC_037335 NC_037334 NC_037332 NC_037331 NC_037328"
for seq in `echo $chr |tr ' ' '\n' ` ; do
bsub -J odgi -n 4 -R "span[hosts=1] rusage[mem=10GB] select[maxmem>200GB]" -o %J.out -e %J.err -q smp \
"
odgi inject -i /public/home/qtyao/WZ_data_transfer/cactus/bos/work/cattle31.${seq}.1.full.og -t 4 -b ${seq}.1.bed -o ${seq}.1.inj.og -t 10
cut -f4 ${seq}.1.bed > ${seq}.1.name
odgi untangle -R ${seq}.1.name -i ${seq}.1.inj.og -g -j 0.1 -n 100 -t 4 > ${seq}.1.untangle.tsv
"
done


```


### 2.Add OR's ID and bubble's length in "all_bubble.ORgt"
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

### 3.Extract OR bubble's syntenic regions across 31 cattle assemblies
The strategy is：

1. Looking for syntenic regions of bubbles for each 31 cattle genomes based on MC graph, if there are satisified projections, take the coordinates;
2. If they are bad projections, take minigraph projection as coordinates;
3. Check bubbles's length, which should equal to projection range as above, if the projection range and bubbles length have a difference over 500 bp and difference/length(bubble) >0.1, discard this projection
```
cd /public/home/qtyao/WZ_data_transfer/minigraph/5.bubble_boundary
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
adj=$(awk -v a=$pos '$2==a' ../../1.minigraph/${chr}/${ind}.bed |awk -F'[_:]' '{print $(NF-2)}')
beg=$(awk -v a=$pos '$2==a' ../../1.minigraph/${chr}/${ind}.bed |awk -F'[_:]' '{print $(NF-1)}')
end=$(awk -v a=$pos '$2==a' ../../1.minigraph/${chr}/${ind}.bed |awk -F'[_:]' '{print $NF}')
beg=$(($beg + $adj - 2000))
end=$(($end + $adj + 2000))
nonchr=$(awk -v a=$pos '$2==a' ../../1.minigraph/${chr}/${ind}.bed |awk -F':' '{print $4}' |sed 's/_[^_]*$//')
strand=$(awk -v a=$pos '$2==a' ../../1.minigraph/${chr}/${ind}.bed |awk -F':' '{print $3}')
gvd=$(grep $ind gtass.tmp |awk '{print "GT"$2}')
orid=$(awk -v a=$vcfid -v b=$gvd '$1==a && $2==b {print $4}' all_bubble.ORgt.withID.withLen)
if [ "$strand" == "-" ]; then
strand=0
elif [ "$strand" == "+" ]; then
strand=1
fi
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
adj=$(awk -v a=$pos '$2==a' ../../1.minigraph/${chr}/${ind}.bed |awk -F'[_:]' '{print $(NF-2)}')
beg=$(awk -v a=$pos '$2==a' ../../1.minigraph/${chr}/${ind}.bed |awk -F'[_:]' '{print $(NF-1)}')
end=$(awk -v a=$pos '$2==a' ../../1.minigraph/${chr}/${ind}.bed |awk -F'[_:]' '{print $NF}')
beg=$(($beg + $adj - 2000))
end=$(($end + $adj + 2000))
nonchr=$(awk -v a=$pos '$2==a' ../../1.minigraph/${chr}/${ind}.bed |awk -F':' '{print $4}' |sed 's/_[^_]*$//')
strand=$(awk -v a=$pos '$2==a' ../../1.minigraph/${chr}/${ind}.bed |awk -F':' '{print $3}')
gvd=$(grep $ind gtass.tmp |awk '{print "GT"$2}')
orid=$(awk -v a=$vcfid -v b=$gvd '$1==a && $2==b {print $4}' all_bubble.ORgt.withID.withLen)
if [ "$strand" == "-" ]; then
strand=0
elif [ "$strand" == "+" ]; then
strand=1
fi
aa=$(awk -v a=$gvd -v b=$vcfid '$1==b && $2==a {print $5}' all_bubble.ORgt.withID.withLen)
echo -e "$nonchr\t$beg\t$end\t$ind\t$strand\t$orid\t$gvd\t$aa\tminigraph"
done |awk 'BEGIN{OFS="\t"} {print ($8 > ($3-$2-4000)) ? ($8 - ($3-$2-4000) ) : (($3-$2-4000) - $8), $0}' |awk '$1<500 || $1/$9 <0.1'  |cut  -f2- >> projection.sort.filter.bed
rm problem gtass.tmp
sort -k6,6 projection.sort.filter.bed > projection.sort.filter.bed0

echo "" > boundary_5000_sequence/${vcfid}.sequence.fasta
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
done  >> boundary_5000_sequence/${vcfid}.sequence.fasta
rm  projection.sort.bed projection.sort.filter.bed0
fi
done
```

### 4. Extract break point ± 100 bp sequence, including two insertion boundaries and one deletion boundary
```
for vcfid in `ls *.mafft.aln |awk -F".sequence" '{print $1}'` ; do
mkdir ${vcfid}
mv ${vcfid}.sequence.fasta ${vcfid}.sequence.fasta.mafft.aln ${vcfid}
cp process_seq.sh boundary_detection.sh run_boundary.sh ${vcfid}
cd /public/home/qtyao/WZ_data_transfer/minigraph/5.bubble_boundary/boundary_2000_sequence/${vcfid}
sed 's/-/N/g' ${vcfid}.sequence.fasta.mafft.aln > ${vcfid}.sequence.fasta.mafft.aln.fasta
~/WZ_data_transfer/tools/minimap2/minimap2 -cx asm20 ../all.coding.cdhit98.rename.fasta ${vcfid}.sequence.fasta.mafft.aln.fasta  > ${vcfid}.sequence.fasta.mafft.aln.fasta.paf 
seqkit seq  ${vcfid}.sequence.fasta.mafft.aln -w 1000000 > ${vcfid}.sequence.fasta.mafft.aln.fasta.linear
awk 'BEGIN{OFS="\t"} {print $1,$3,$4}' ${vcfid}.sequence.fasta.mafft.aln.fasta.paf  |bedtools sort -i - |bedtools merge -i - |\
awk -v a=$vcfid 'BEGIN{OFS="\t"} {print a,$2,$3,$1}' |bedtools sort -i - |bedtools merge -i - -c 4 -o collapse | awk 'BEGIN{OFS="\t"} $3-$2 > 700 {print $1,$2,$3,$4}' > ${vcfid}.OR.bed
grep "^>" ${vcfid}.sequence.fasta |awk -F'>' '{print $2}' > ${vcfid}.sequence.fasta.name
cat ${vcfid}.OR.bed |while IFS= read -r line; do
del=$(echo $line |awk '{print $4}' |tr ',' '\n' | grep -v -f - ${vcfid}.sequence.fasta.name |tr '\n' ',' |sed 's/,$//')
echo $line $del
done |awk '$4 !="" && $5 !=""' > ${vcfid}.OR.bed.anno
bsub -J boundary -n 2 -R "span[hosts=1] rusage[mem=4GB]" -o %J.out -e %J.err -q normal \
"
sh ./run_boundary.sh ${vcfid}
"
cd /public/home/qtyao/WZ_data_transfer/minigraph/5.bubble_boundary/boundary_2000_sequence
done
```

### 5. Calculate boundary's simialarity
```
for id in `cut -f1 combined_boundary.seq.fai |awk -F':' '{print $2":"$3":"$4}' |uniq`; do

ins1up=$(grep ${id} combined_boundary.seq.fai |grep "^ins_b1" |cut -f1 | ~/WZ_data_transfer/tools/seqtk-master/seqtk subseq combined_boundary.seq - |tail -n1 |cut -c 51-100)
ins2up=$(grep ${id} combined_boundary.seq.fai |grep "^ins_b2" |cut -f1 | ~/WZ_data_transfer/tools/seqtk-master/seqtk subseq combined_boundary.seq - |tail -n1 |cut -c 51-100)
delup=$(grep ${id} combined_boundary.seq.fai |grep "^del" |cut -f1 | ~/WZ_data_transfer/tools/seqtk-master/seqtk subseq combined_boundary.seq -    |tail -n1 |cut -c 51-100)

ins1dw=$(grep ${id} combined_boundary.seq.fai |grep "^ins_b1" |cut -f1 | ~/WZ_data_transfer/tools/seqtk-master/seqtk subseq combined_boundary.seq - |tail -n1 |cut -c 101-150)
ins2dw=$(grep ${id} combined_boundary.seq.fai |grep "^ins_b2" |cut -f1 | ~/WZ_data_transfer/tools/seqtk-master/seqtk subseq combined_boundary.seq - |tail -n1 |cut -c 101-150)
deldw=$(grep ${id} combined_boundary.seq.fai |grep "^del" |cut -f1 | ~/WZ_data_transfer/tools/seqtk-master/seqtk subseq combined_boundary.seq -    |tail -n1 |cut -c 101-150)

ins1upn=$(echo $ins1up |grep -o "n" |wc -l)
ins2upn=$(echo $ins2up |grep -o "n" |wc -l)
delupn=$(echo $delup |grep -o "n" |wc -l )
ins1dwn=$(echo $ins1dw |grep -o "n" |wc -l)
ins2dwn=$(echo $ins2dw |grep -o "n" |wc -l)
deldwn=$(echo $deldw |grep -o "n" |wc -l )

if [ ${#ins1up} -gt 45 ] && [ ${#ins2up} -gt 45 ] && [ ${#delup} -gt 45 ] && [ ${#ins1dw} -gt 45 ] && [ ${#ins2dw} -gt 45 ] && [ ${#deldw} -gt 45 ] && [ ${#ins1upn} -lt 5 ] && [ ${#ins2upn} -lt 5 ] && [ ${#delupn} -lt 5 ] && [ ${#ins1dwn} -lt 5 ] && [ ${#ins2dwn} -lt 5 ] && [ ${#deldwn} -lt 5 ]; then
echo -e ">ins1up\n$ins1up\n>delup\n$delup" > ins1up_delup.fasta
mafft ins1up_delup.fasta  > ins1up_delup.fasta.aln
match=$(seqkit seq  ins1up_delup.fasta.aln  -w 999 |grep -v "^>" | awk 'NR==1 {line1=$0} NR==2 {line2=$0} 
     END {
         count=0
         len=length(line1)
         for (i=1; i<=len; i++) {
             if (substr(line1, i, 1) == substr(line2, i, 1)) {
                 count++
             }
         }
         print count
     }' -)
length=$(seqkit seq  ins1up_delup.fasta.aln  -w 999 |grep -v "^>" |head -1 |awk '{print length($0)}' )  
ins1up_delupsim=$(echo "scale=2; $match / $length" | bc)

echo -e ">ins1dw\n$ins1dw\n>deldw\n$deldw" > ins1dw_deldw.fasta
mafft ins1dw_deldw.fasta  > ins1dw_deldw.fasta.aln
match=$(seqkit seq  ins1dw_deldw.fasta.aln  -w 999 |grep -v "^>" | awk 'NR==1 {line1=$0} NR==2 {line2=$0} 
     END {
         count=0
         len=length(line1)
         for (i=1; i<=len; i++) {
             if (substr(line1, i, 1) == substr(line2, i, 1)) {
                 count++
             }
         }
         print count
     }' -)
length=$(seqkit seq  ins1dw_deldw.fasta.aln  -w 999 |grep -v "^>" |head -1 |awk '{print length($0)}' )  
ins1dw_deldwsim=$(echo "scale=2; $match / $length" | bc)

echo -e ">ins2up\n$ins2up\n>delup\n$delup" > ins2up_delup.fasta
mafft ins2up_delup.fasta  > ins2up_delup.fasta.aln
match=$(seqkit seq  ins2up_delup.fasta.aln  -w 999 |grep -v "^>" | awk 'NR==1 {line1=$0} NR==2 {line2=$0} 
     END {
         count=0
         len=length(line1)
         for (i=1; i<=len; i++) {
             if (substr(line1, i, 1) == substr(line2, i, 1)) {
                 count++
             }
         }
         print count
     }' -)
length=$(seqkit seq  ins2up_delup.fasta.aln  -w 999 |grep -v "^>" |head -1 |awk '{print length($0)}' )  
ins2up_delupsim=$(echo "scale=2; $match / $length" | bc)


echo -e ">ins2dw\n$ins2dw\n>deldw\n$deldw" > ins2dw_deldw.fasta
mafft ins2dw_deldw.fasta  > ins2dw_deldw.fasta.aln
match=$(seqkit seq  ins2dw_deldw.fasta.aln  -w 999 |grep -v "^>" | awk 'NR==1 {line1=$0} NR==2 {line2=$0} 
     END {
         count=0
         len=length(line1)
         for (i=1; i<=len; i++) {
             if (substr(line1, i, 1) == substr(line2, i, 1)) {
                 count++
             }
         }
         print count
     }' -)
length=$(seqkit seq  ins2dw_deldw.fasta.aln  -w 999 |grep -v "^>" |head -1 |awk '{print length($0)}' )  
ins2dw_deldwsim=$(echo "scale=2; $match / $length" | bc)

echo $id $ins1up_delupsim $ins1dw_deldwsim $ins2up_delupsim $ins2dw_deldwsim
else
  echo -e "$id : Bad boundary"
fi
done > combined_boundary_similarity_50bp.results
```
