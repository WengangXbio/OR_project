fastapath=$1
OR_path=$2
CHROM=$3
genomelist=$(cut -f1 ${fastapath}/${CHROM}/${CHROM}.pggb.fasta.fai |awk -F'#' '{print $1}' |sort |uniq)
for genome in `echo $genomelist`; do
OR_dir=$(ls -l $OR_path | grep '^d' | awk '{print $9}' |grep ${genome})
awk -F'_' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' ${OR_path}/${OR_dir}/ORannotation_itera2_final_dead.fasta.fai |\
sort -k1,1 -k2,2n > pseudo.bed
sh modify_coor.sh pseudo.bed modify.pseudo.bed
sed "s/@#@/_/g" modify.pseudo.bed > pseudo.bed
cut -f1 ${fastapath}/${CHROM}/${CHROM}.pggb.fasta.fai | grep ${genome} > seqfile
./precess_bed.sh seqfile pseudo.bed
done > dead_annotation.bed
rm  pseudo.bed modify.pseudo.bed seqfile

