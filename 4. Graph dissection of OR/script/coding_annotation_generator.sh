fastapath=$1
OR_path=$2
CHROM=$3
genomelist=$(cut -f1 ${fastapath}/${CHROM}/${CHROM}.pggb.fasta.fai |awk -F'#' '{print $1}' |sort |uniq)
for genome in `echo $genomelist`; do
OR_dir=$(ls -l $OR_path | grep '^d' | awk '{print $9}' |grep ${genome})
awk -F'_' 'BEGIN{OFS="\t"} {print $2,$3,$4,$5}' ${OR_path}/${OR_dir}/ORannotation_itera2_final_func_dna_ORs.fasta.fai |\
sort -k1,1 -k2,2n > coding.bed
sh modify_coor.sh coding.bed modify.coding.bed
sed "s/@#@/_/g" modify.coding.bed > coding.bed
cut -f1 ${fastapath}/${CHROM}/${CHROM}.pggb.fasta.fai | grep ${genome} > seqfile
./precess_bed.sh seqfile coding.bed
done > coding_annotation.bed
rm coding.bed modify.coding.bed seqfile
