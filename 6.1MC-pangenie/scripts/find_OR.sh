suffix=$1
cat cattle31.raw.vcfbub.vcf.diploid.noheader.removeN.vcf.${suffix} |while IFS= read -r line; do
vcfid=$(echo $line |awk '{print $1"_"$2}')
echo $line |awk '{print $4}' > seq.fasta
echo $line |awk '{print $5}'|tr ',' '\n' >> seq.fasta
awk '{print ">GT"NR-1,$0}' seq.fasta |tr ' ' '\n' >  GTseq.fasta
samtools faidx GTseq.fasta
maxl=$(cut -f2 GTseq.fasta.fai |sort |tail -1)
if [ $maxl -gt 500 ]; then
minimap2 -cx asm5 all.coding.cdhit98.rename.fasta GTseq.fasta  > GTseq.fasta.paf
pafl=$(cat GTseq.fasta.paf  |wc -l)
if [ "$pafl" != 0 ]; then
awk 'BEGIN{OFS="\t"} {print $1,$3,$4,$6}' GTseq.fasta.paf  |sort -k1,1 -k2,2n -k4,4 |bedtools merge -i - -c 4 -o collapse -delim "|" |awk '$3-$2>700' > aln.gt
tgl=$(cat aln.gt |wc -l)
if [ "$tgl" != 0 ]; then
cat aln.gt |while IFS= read -r gtl; do
gta=$(echo $gtl |awk '{print $1}' )
sss=$(echo $gtl |awk '{print $2}' )
eee=$(echo $gtl |awk '{print $3}' )
ooo=$(echo $gtl |awk '{print $4}' )
mm=$(echo $ooo |tr '|' '\n' |sort |tr '\n' '|')
echo -e "$gta\t$sss\t$eee\t$mm"
done > aln.gt0
grep '^>' GTseq.fasta |awk -F'>' '{print $2}' | while IFS= read -r gt; do
vb=$(awk -v a=$gt '$1==a'  aln.gt0)
if [ -n "$vb" ]; then
awk -v a=$gt '$1==a' aln.gt0 |cut -f4 |sed ':a;N;$!ba;s/\n/),(/g' |awk -v a=$gt -v b=$vcfid 'BEGIN{OFS="\t"} {print b,a,"("$0")"}'
else
echo -e "$vcfid\t$gt\t*"
fi
done
rm aln.gt0
fi
rm aln.gt
fi
rm GTseq.fasta.paf
fi
rm seq.fasta GTseq.fasta GTseq.fasta.fai
done > ${suffix}.bubble.ORgt
