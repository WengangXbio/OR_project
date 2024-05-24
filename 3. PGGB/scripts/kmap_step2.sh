#Extract chromosomes/scaffolds/cotigs from assemblies for chromosome-pangenome
wd="/share/home/yzwl_zhangwg/OR_project/PGGB/pan_by_chr_new"
for i in `seq 1 31`  ; do
chrn=chr$i
mkdir $wd/$chrn -p
ass=$(awk -v a=$chrn '$1==a {print $2}' ARS20_chr_NC.info)
alignmap=$(ls *.align_map |grep -v GCA_007844835.1_NRC_Mithun_1_genomic.fna |grep -v GCF_000754665.1_Bison_UMD1.0_genomic.fna)

for fna in `echo $alignmap`; do
echo -e "processing $chrn $fna"
grep $ass $fna | while IFS= read -r line; do
fa=$(echo $line |awk '{print $1}')
cg=$(echo $line |awk '{print $2}')
st=$(echo $line |awk '{print $3}')
ed=$(echo $line |awk '{print $4}')
fasta="/share/home/yzwl_zhangwg/OR_project/PGGB/assembly/${fa}"
samplename=$(echo $fa|awk -F'.' '{print $1}')
if [ "$ed" = 0 ]
then
samtools faidx ${fasta} ${cg} > $wd/$chrn/${cg}_${st}.fasta
newnames=">${samplename}#1#${cg}"
sed -i "1s/.*/$newnames/" $wd/$chrn/${cg}_${st}.fasta 
else
echo -e "$cg\t$st\t$ed" |bedtools getfasta -fi ${fasta} -bed - > $wd/$chrn/${cg}_${st}.fasta
newnames=">${samplename}#1#${cg}_${st}_${ed}"
sed -i "1s/.*/$newnames/" $wd/$chrn/${cg}_${st}.fasta 
fi
done
done
cat $wd/$chrn/*.fasta > $wd/$chrn/${chrn}.pggb.fasta
samtools faidx $wd/$chrn/${chrn}.pggb.fasta
done
