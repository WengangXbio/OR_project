chrn=$1
ass=$2
fnalist=$(ls *.paf |awk -F'.paf' '{print $1}')
wd="/share/home/yzwl_zhangwg/OR_project/PGGB/pan_by_chr_new_new"
for fna in `echo $fnalist`; do
echo -e "processing $chrn $fna"
grep $ass ${fna}/${fna}.align_map | while IFS= read -r line; do
fa=$(echo $line |awk '{print $1}')
cg=$(echo $line |awk '{print $2}')
st=$(echo $line |awk '{print $3}')
ed=$(echo $line |awk '{print $4}')
fasta="/share/home/yzwl_zhangwg/OR_project/PGGB/assembly/${fa}"
samplename=$(echo $fa|awk -F'.' '{print $1}')
if [ "$ed" = 0 ]
then
samtools faidx ${fasta} ${cg} > $wd/$chrn/${fna}_${cg}_${st}.fasta
newnames=">${samplename}#1#${cg}"
sed -i "1s/.*/$newnames/" $wd/$chrn/${fna}_${cg}_${st}.fasta 
else
echo -e "$cg\t$st\t$ed" |bedtools getfasta -fi ${fasta} -bed - > $wd/$chrn/${fna}_${cg}_${st}.fasta
newnames=">${samplename}#1#${cg}_${st}_${ed}"
sed -i "1s/.*/$newnames/" $wd/$chrn/${fna}_${cg}_${st}.fasta
fi
done
done
#cat $wd/$chrn/*.fasta > $wd/$chrn/${chrn}.pggb.fasta
#samtools faidx $wd/$chrn/${chrn}.pggb.fasta


