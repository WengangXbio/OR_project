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
```
### Implement minigraph by chromosome
```
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

### Call SV on each assembly
```
for ctg in `seq 37328 37357 ` ; do
cd  NC_0${ctg}.1
chr=NC_0${ctg}.1
gfa=${chr}.gfa
for genome in `ls /public/home/qtyao/WZ_data_transfer/minigraph/0.assembly/ |grep fasta |grep -v fai |awk -F'.' '{print $1}'` ; do
bsub -J bed -n 12 -R "span[hosts=1] rusage[mem=2GB] " -o %J.out -e %J.err -q normal \
"
/public/home/qtyao/WZ_data_transfer/tools/minigraph-master/minigraph -cxasm --call -t12 ${gfa} ${genome}_NC_0${ctg}.1.fasta > ${genome}.bed
"
done
cd /public/home/qtyao/WZ_data_transfer/minigraph/1.minigraph
done
```
### Merge SV for all assemblies to generate minigraph-based VCF
```
for ctg in `seq 37328 37357 ` ; do
cd  NC_0${ctg}.1
chr=NC_0${ctg}.1
awk '{print $0".bed"}' ../genome_prefix |tr '\n' ' ' |xargs paste | /public/home/qtyao/WZ_data_transfer/tools/mg-cookbook-v1_x64-linux/k8 \
/public/home/qtyao/WZ_data_transfer/tools/mg-cookbook-v1_x64-linux/mgutils.js merge -s <(cat ../genome_prefix ) - | gzip > NC_0${ctg}.1.sv.bed.gz
/public/home/qtyao/WZ_data_transfer/tools/mg-cookbook-v1_x64-linux/k8 /public/home/qtyao/WZ_data_transfer/tools/mg-cookbook-v1_x64-linux/mgutils-es6.js \
merge2vcf -r0 NC_0${ctg}.1.sv.bed.gz > NC_0${ctg}.1.sv.vcf
grep "^##" NC_0${ctg}.1.sv.vcf | cat - ../replace_title - <(grep -v "^#" NC_0${ctg}.1.sv.vcf) > NC_0${ctg}.1.sv.vcf0
cd /public/home/qtyao/WZ_data_transfer/minigraph/1.minigraph
done
```

### Modifing vcf for each chromosome
```
for ctg in `seq 37328 37357 ` ; do
chr=NC_0${ctg}.1
mkdir ${chr}
cd ${chr}
ln -s /public/home/qtyao/WZ_data_transfer/minigraph/1.minigraph/${chr}/${chr}.sv.vcf0 ${chr}.sv.vcf0
ln -s /public/home/qtyao/WZ_data_transfer/minigraph/1.minigraph/${chr}/${chr}.gfa ${chr}.gfa 
ln -s /public/home/qtyao/WZ_data_transfer/minigraph/2.vgvcf/vcf_generation.sh vcf_generation.sh
grep '^S' ${chr}.gfa > out.gfa.node
bsub -J bed -n 2 -R "span[hosts=1] rusage[mem=4GB] " -o %J.out -e %J.err -q normal \
"
sh ./vcf_generation.sh ${chr}.sv.vcf0 out.gfa.node GCF_002263795  # input_vcf, input_node, reference_id
"
cd /public/home/qtyao/WZ_data_transfer/minigraph/2.vgvcf
done
```

### Correcting chromosome name of modified vcf
```
for ctg in `seq 37328 37357 ` ; do
chr=NC_0${ctg}.1
chr0=NC_0${ctg}.1_0
cd ${chr}
sed "s/${chr0}/${chr}/g" vcf.body >  ${chr}.sv.vcf0.vg.vcf0
cut -f1 ${chr}.sv.vcf0.vg.vcf0 |uniq
cd /public/home/qtyao/WZ_data_transfer/minigraph/2.vgvcf
done
```

### Find inconsistent site (including inconsistent sites and duplicated node)
Inconsistent sites: The modified vcf's reference seqence is not consistent with reference sequence

Duplicated node: Some of sequences of two alternative genotypes are exactly same  
```
for ctg in `seq 37328 37357` ; do
chr=NC_0${ctg}.1
cd ${chr}
echo > conflict_genotype
echo > inconsistent_with_ref
cp ../vcf_correction.sh ./
bsub -J correction -n 2 -R "span[hosts=1] rusage[mem=4GB] " -o %J.out -e %J.err -q normal \
"
sh ./vcf_correction.sh ${chr}
"
cd /public/home/qtyao/WZ_data_transfer/minigraph/2.vgvcf
done
```
