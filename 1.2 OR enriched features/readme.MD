## Genome-wide profiling based on pangenome and featured with OR regions


### window depth by 10K window project on ARS2.0 
```
for nchr in `seq 1 29`  ; do
odgi build -t 4 -P -g chr${nchr}.gfa -o chr${nchr}.og -O
rhe=$(echo chr${nchr}) 
che=$(awk -v a=${rhe} '$1==a {print $2}' chr_anno.txt)
odgi depth -i ${rhe}.og -r ${che}#1#${rhe} | bedtools makewindows -b /dev/stdin -w 10000 > window_depth/${rhe}.w10k.bed
odgi depth -i ${rhe}.og -b window_depth/${rhe}.w10k.bed --threads 2 | bedtools sort > window_depth/${rhe}.w10k.depth.bed
done
```

### variants density on CDS regions by 10K window walking
```
awk -F'\t' 'BEGIN{OFS="\t"} $3=="CDS" {print $1,$4,$5}' GCF_002263795.3_ARS-UCD2.0_genomic.gff > GCF_002263795.3_ARS-UCD2.0_genomic.gff.CDS.bed
for nchr in `seq 1 29`  ; do
chr="chr${nchr}"
ass=$(awk -v a=${chr} '$1==a {print $2}' chr_anno.txt)
sed -i "s/${ass}/${chr}/g" GCF_002263795.3_ARS-UCD2.0_genomic.gff.CDS.bed
done
grep "chr" GCF_002263795.3_ARS-UCD2.0_genomic.gff.CDS.bed > GCF_002263795.3_ARS-UCD2.0_genomic.gff.CDS.bed0
for nchr in `seq 1 29`  ; do
chr="chr${nchr}"
length=$(awk -v a=${chr} '$1==a {print $3}' chr_anno.txt)
grep -v "^#" /share/home/yzwl_zhangwg/OR_project/PGGB/work/${chr}.vcf | awk -v a=${chr} 'BEGIN{OFS="\t"} {print a,$2,$2+1}' - > tmp.vcf.bed
bedtools intersect -b tmp.vcf.bed -a GCF_002263795.3_ARS-UCD2.0_genomic.gff.CDS.bed0 > tmp.CDS.vcf.bed
echo -e "$chr\t"0"\t$length" |bedtools makewindows -b - -w 1000 > tmp.window.bed
bedtools intersect -b tmp.CDS.vcf.bed -a tmp.window.bed -c > window_vcf_density/${chr}.CDS.vcf.density
rm tmp.vcf.bed tmp.CDS.vcf.bed tmp.window.bed
done
```

### variants density by 10K window walking
```
for nchr in `seq 1 29`  ; do
chr="chr${nchr}"
length=$(awk -v a=${chr} '$1==a {print $3}' chr_anno.txt)
grep -v "^#" /share/home/yzwl_zhangwg/OR_project/PGGB/work/${chr}.vcf | awk -v a=${chr} 'BEGIN{OFS="\t"} {print a,$2,$2+1}' - > tmp.vcf.bed
echo -e "$chr\t"0"\t$length" |bedtools makewindows -b - -w 10000 > tmp.window.bed
bedtools intersect -b tmp.vcf.bed -a tmp.window.bed -c > window_vcf_density/${chr}.vcf.w10k.density
rm tmp.window.bed tmp.vcf.bed
done
```

### PAV average value by 10K window walking
```
cp /share/home/yzwl_zhangwg/OR_project/PGGB/work/window_depth/*.w10k.bed ./window_pav
cd /share/home/yzwl_zhangwg/OR_project/PGGB/work/
for nchr in `seq 1 29`  ; do
chr="chr${nchr}"
odgi pav -i ${chr}.og -b window_pav/${chr}.w10k.bed > window_pav/${chr}.w10k.pav
done
```

### Plot in R
```
library("karyoploteR")
library(rtracklayer)
setwd("E:\\OR_project\\GCF_002263795.3_ARS-UCD2.0_genomic.gff")
karyotype=read.table("kayotype.txt",head=F)
gff.file="GCF_002263795.3_ARS-UCD2.0_genomic.OR.gff"
features <- import(gff.file)
BS_genome <- toGRanges(karyotype)
kp <- plotKaryotype(genome=BS_genome)
table(features$type)
genes <- features[features$type=="gene" ]
pseudogene <- features[features$type=="pseudogene"]
pp <- getDefaultPlotParams(plot.type=2)
pp$data1outmargin <- 100
pp$data2outmargin <- 100
pp$topmargin <- 450
kp <- plotKaryotype(genome=BS_genome, ideogram.plotter = NULL, plot.type=2, plot.params = pp)
kpAddCytobandsAsLine(kp)
kpAddMainTitle(kp, "Olfactory receptor genes or pseudogenes", cex=2)
kpPlotRegions(kp, data=pseudogene, avoid.overlapping = FALSE, col="gold")
kpPlotRegions(kp, data=genes, avoid.overlapping = FALSE, col="deepskyblue")

for(i in 1:29){
#depth=read.table(paste("chr",i,".w10k.depth.bed",sep=""),head=F,comment.char="")
#x=depth$V3
#y=depth$V4-median(depth$V4)
#y0=y[-which(y<=0)]+1
#x0=x[-which(y<=0)]
#y0=log2(y0)/max(log2(y0))
#kpPoints(kp, chr=paste("NC_0373",i+27,".1",sep=""), x=x0, y=y0, data.panel=2,cex=0.1)
#vcfden=read.table(paste("chr",i,".vcf.density",sep=""),head=F,comment.char="")
#x=vcfden$V3
#y=vcfden$V4
#y0=y/max(y)
#x0=x
#kpPoints(kp, chr=paste("NC_0373",i+27,".1",sep=""), x=x0, y=y0, data.panel=2,cex=0.1)

#pav=read.table(paste("chr",i,".w10k.pav",sep=""),head=T,comment.char="")
#x=c()
#y=c()
#for(j in unique(pav$end)){
#x=c(x,j)
#y=c(y,mean(pav[which(pav$end==j),]$pav))
#}
#y0=1-y
#x0=x
#kpPoints(kp, chr=paste("NC_0373",i+27,".1",sep=""), x=x0, y=y0, data.panel=2,cex=0.1)

vcfden=read.table(paste("chr",i,".vcf.w10k.density",sep=""),head=F,comment.char="")
x=vcfden$V3
y=vcfden$V4
y0=y/max(y)
x0=x
kpPoints(kp, chr=paste("NC_0373",i+27,".1",sep=""), x=x0, y=y0, data.panel=2,cex=0.1)
}


```
