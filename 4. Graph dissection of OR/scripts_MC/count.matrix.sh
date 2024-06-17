nchr=$1
NC=$2
cp ../genomelist ./
itera=1
cp ${NC}.func.proc.bed itera.coding_annotation.bed
cp ${NC}.dead.proc.bed dead_annotation.bed
lw=$(awk 'BEGIN{OFS="\t"} {print $1,$3,$4,$2}' injref.untangle.tsv.collapse | bedtools intersect -a itera.coding_annotation.bed -b - -wa -wb |wc -l)
genome=GCF_002263795
echo > target.gene
while [ $lw -ne 0 ] ; do
grep ${genome} itera.coding_annotation.bed |cut -f4 > itera.ORlist
./pseudo_matrix.sh itera.ORlist genomelist injref.untangle.tsv.collapse itera.coding_annotation.bed dead_annotation.bed ${itera}
./OR_counter.sh itera.ORlist genomelist injref.untangle.tsv.collapse itera.coding_annotation.bed ${itera}
itera=$((itera + 1))
grep -v -f target.gene itera.coding_annotation.bed > temp 
mv temp itera.coding_annotation.bed
genome=$(awk 'BEGIN{OFS="\t"} {print $1,$3,$4,$2}' injref.untangle.tsv.collapse | bedtools intersect -a itera.coding_annotation.bed -b - -wa -wb |\
cut -f1-4| awk -F'#' '{print $1}' |sort |uniq -c |sort -k1,1nr |awk '{print $2}' |head -n1)
lw=$(awk 'BEGIN{OFS="\t"} {print $1,$3,$4,$2}' injref.untangle.tsv.collapse | bedtools intersect -a itera.coding_annotation.bed -b - -wa -wb |wc -l)
#echo $itera
done
rm  itera.* target.gene intersect.dead test.intersect.dead
cat or.counts_*.matrix | awk -v chr=${nchr} '{print chr"."$0}' > or.counts.matrix
cat or.pseudo_*.matrix | awk -v chr=${nchr} '{print chr"."$0}' > or.pseudo.matrix
sed '1i\\' genomelist | tr '\n' '\t' > header
printf '\n' >> header
cat header or.counts.matrix > or.counts.matrix0
cat header or.pseudo.matrix > or.pseudo.matrix0

