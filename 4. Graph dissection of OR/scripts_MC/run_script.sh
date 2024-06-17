####split og file by chromosome
bsub -J repeatmask -n 4 -R "span[hosts=1] rusage[mem=4GB] select[maxmem>224800]" -o %J.out -e %J.err -q normal '
sh split_chr.sh
'
####writeout path names
bsub -J repeatmask -n 4 -R "span[hosts=1] rusage[mem=4GB] select[maxmem>224800]" -o %J.out -e %J.err -q normal '
sh write_path.sh
'

####modify annotation bed file
cd or_anno_transfer/

for genome in `ls *.g2or_itera2_func.bed |awk -F'.' '{print $1}' |sort |uniq`; do
func=$(ls *g2or_itera2_func.bed |grep ${genome} |grep -v "dead")
dead=$(ls *g2or_itera2_func.bed |grep ${genome} |grep "dead")
awk -F'\t' -v a=${genome} 'BEGIN{OFS="\t"}{print a"#0#"$1"#0",$2,$3,a"_"$1"_"$2}' ${func} > ${genome}.func.bed
awk -F'\t' -v a=${genome} 'BEGIN{OFS="\t"}{print a"#0#"$1"#0",$2,$3,a"_"$1"_"$2}' ${dead} > ${genome}.dead.bed
done




####process bed file by adjusting coordinates
bsub -J procbed -n 4 -R "span[hosts=1] rusage[mem=6GB] select[maxmem>100GB]" -o %J.out -e %J.err -q normal '
sh procbed.sh
'

for NC in `cut -f1 ref.info |awk -F"#" '{print $3}'` ; do
if [ -s "./${NC}.func.proc.bed" ]; then
mkdir chr_${NC} -p
cd chr_${NC}
ln -s /public/home/qtyao/WZ_data_transfer/cactus/bos/work/cattle31.${NC}.full.og ./
cp /public/home/qtyao/WZ_data_transfer/cactus/bos/work/${NC}.dead.proc.bed ./
cp /public/home/qtyao/WZ_data_transfer/cactus/bos/work/${NC}.func.proc.bed ./
cp /public/home/qtyao/WZ_data_transfer/cactus/bos/work/projection.sh ./
sed -i "s/OMG/${NC}/g" projection.sh
bsub < projection.sh
cd ..
else 
   echo "${NC} is blank"
fi 
done


