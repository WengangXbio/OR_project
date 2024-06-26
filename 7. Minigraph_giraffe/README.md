## Minigraph + giraffe 
Using minigraph build graph and then mapping short reads to call SV
### Seperate sequence by chromosome
```
cd /public/home/qtyao/WZ_data_transfer/minigraph/1.minigraph
for ctg in `seq 37328 37357 |awk '{print} END {print 82638}'` ; do
mkdir NC_0${ctg}.1
cd  NC_0${ctg}.1
for genome in `ls /public/home/qtyao/WZ_data_transfer/minigraph/0.assembly/ |grep fasta |grep -v fai |awk -F'.' '{print $1}'` ; do
if [[ ${genome} != "GCF_002263795" ]]; then
awk -F'[#:-]' -v a=${genome} 'BEGIN{OFS="\t"} $1==a {print $3,$5,$6,$3"_"$5}' ../cattle31.NC_0${ctg}.1.syntenic.regions | \
bedtools getfasta -fi /public/home/qtyao/WZ_data_transfer/minigraph/0.assembly/${genome}.fasta  -bed - -nameOnly > ${genome}_NC_0${ctg}.1.fasta
fi
if [[ ${genome} == "GCF_002263795" ]]; then
awk -F'[#:-]' -v a=${genome} 'BEGIN{OFS="\t"} $1==a {print $3,$4,$5,$3"_"$4}' ../cattle31.NC_0${ctg}.1.syntenic.regions | \
bedtools getfasta -fi /public/home/qtyao/WZ_data_transfer/minigraph/0.assembly/${genome}.fasta  -bed - -nameOnly > ${genome}_NC_0${ctg}.1.fasta
fi
done
cd /public/home/qtyao/WZ_data_transfer/minigraph/1.minigraph
done
#implement minigraph by chromosome
for ctg in `seq 37328 37357 |awk '{print} END {print 82638}'` ; do
cd  NC_0${ctg}.1
chr=NC_0${ctg}.1
awk -v a=${chr} '{print $0"_"a".fasta"}' ../genome_prefix |tr '\n' ' ' |\
awk -v a=$chr '{print "/public/home/qtyao/WZ_data_transfer/tools/minigraph-master/minigraph -cxggs -t24 "$0" > "a".gfa" }' > minigraph.sh
bsub -J minigraph -n 24 -R "span[hosts=1] rusage[mem=6GB] " -o %J.out -e %J.err -q normal \
'
sh minigraph.sh
'
cd /public/home/qtyao/WZ_data_transfer/minigraph/1.minigraph
done
```
