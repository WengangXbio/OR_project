chr=$1
cat ${chr}.sv.vcf0.vg.vcf0 |while IFS= read -r line; do
vref=$(echo $line |awk '{print $4}')
valt=$(echo $line |awk '{print $5}')
chr=$(echo $line |awk '{print $1}')
pos=$(echo $line |awk '{print $2}')
name=$(echo $line |awk '{print $3}')
length=${#vref}
end=$((pos + length -1))
gt=$(echo $line |cut -d' ' -f10-)
rref=$(samtools faidx ../../0.assembly/GCF_002263795.fasta ${chr}:${pos}-${end} | grep -v "^>" |tr '[:lower:]' '[:upper:]' |tr -d '\n' )

if [ $vref != $rref ]; then
echo $chr $pos >> inconsistent_with_ref
#echo $line | awk -v a=${rref} -F' ' 'BEGIN{OFS="\t"} $4=a {print $0}'
else

echo $vref,$valt |tr ',' '\n' |awk '{print NR-1,$0}' > pad
awk '{print $1}' pad |while IFS= read -r gtng; do
ind=$(echo $gt |tr ' ' '\n' |paste - ../genome_prefix |awk -v a=$gtng '$1==a {print $2}' |tr '\n' ',' |sed 's/,$//')
awk -v a=$gtng -v b=$ind '{ if ($1 == a) {print $0,b}else{print $0}}' pad > temp
mv temp pad
done
sequence=$(awk '{print $2}' pad |sort |uniq -c |awk '$1==2 {print $2}' |head -1)
if [ -n "$sequence" ]; then
while [ -n "$sequence" ] ; do
a1=$(awk -v a=$sequence '$2==a {print $1}' pad | head -1)
awk -v a=$sequence '$2==a {print $1}' pad | sed '1d'|while IFS= read -r a2; do
a2ind=$(awk -v a=$a2 '$1==a {print $3}' pad)
awk -v a=$a1 -v b=$a2ind '{ if ($1 == a) {print $0","b}else{print $0}}' pad > temp
mv temp pad
awk -v a=$a2 '{ if ($1 != a) {print $0}}' pad |awk '{print NR-1,$2,$3}'   > temp
mv temp pad
done
sequence=$(awk '{print $2}' pad |sort |uniq -c |awk '$1>=2 {print $2}' |head -1)
done
sed -i 's/,,/,/g' pad
newgt=$(for ind in `cat ../genome_prefix` ; do
dgs=$(grep $ind pad |awk '{print $1}')
if [ -n "$dgs" ]; then
echo $dgs
else
echo "."
fi
done)
newalt=$(sed '1d' pad |awk '{print $2}'|tr '\n' ','|sed 's/,$//')
echo $chr $pos $name $vref $newalt "60 . AC=1;AF=0;AN=4;NS=4 GT" $newgt |tr ' ' '\t'
echo $chr $pos >> conflict_genotype
else
echo $line | tr ' ' '\t'
fi
fi
done > ${chr}.sv.vcf0.vg.correction.vcf0
