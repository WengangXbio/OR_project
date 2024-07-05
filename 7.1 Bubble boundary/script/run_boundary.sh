vcfid=$1
cat "${vcfid}.OR.bed.anno" | while IFS= read -r line; do
orid=$(echo "$line" | awk '{print $1"_"$2}')
start=$(echo "$line" | awk '{print $2}')
end=$(echo "$line" | awk '{print $3}')
gainid=$(echo "$line" | awk '{print $4}' | tr ',' '\n' | awk -F'_' '{print $(NF-1)}' | sort | uniq)
omitid=$(echo "$line" | awk '{print $5}' | tr ',' '\n' | awk -F'_' '{print $(NF-1)}' | sort | uniq)
gainidlist=$(for id in $gainid; do
echo "$line" | awk '{print $4}' | tr ',' '\n' | grep "$id" | shuf | head -1
done)
omitidlist=$(for id in $omitid; do
echo "$line" | awk '{print $5}' | tr ',' '\n' | grep "$id" | shuf | head -1
done)
for gain in $gainidlist; do
for omit in $omitidlist; do
fnl=$(samtools faidx "${vcfid}.sequence.fasta.mafft.aln.fasta.linear" "${omit}:${start}-${end}" | sed '1d' | tr -d '-' | wc -c)
if (( fnl > 30 )); then
continue
else
sh ./boundary_detection.sh ${vcfid} ${gain} ${omit} ${orid} ${start} ${end}
fi
done
done
done > ${vcfid}.boundary.seq
