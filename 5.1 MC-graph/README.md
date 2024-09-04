## MC PAV or dosage or so forth
cd /home/wgzhang/cactus/bos/new_work

### 1. find orthology gene based on pangeome
```
cd /home/wgzhang/cactus/bos/new_work
cat seq.txt |while IFS= read -r proc; do
echo $proc
chr=$(echo $proc |awk '{print $1}')
nc=$(echo $proc |awk '{print $2}')
awk 'BEGIN{OFS="\t"} {print $1,$3,$4,$2}' /home/wgzhang/cactus/bos/new_work/chr_${nc}/injref.untangle.tsv.collapse |bedtools sort -i - | bedtools merge -i - -c 4 -o collapse |awk '{print $4}'> chr${chr}.group
while true; do
    rep=$(cat chr${chr}.group | tr ',' '\n' |grep -v "^$" | sort | uniq -c | awk '$1 > 1 {print $2}' | head -1)
    if [ -n "$rep" ]; then 
        grep "$rep" chr${chr}.group | tr ',' '\n' | sort | uniq | tr '\n' ',' > rep_add && echo >> rep_add
        grep -v "$rep" chr${chr}.group > rep_removal
        cat rep_add rep_removal > chr${chr}.group
    else
        break
    fi
done
done
```

### 2. count dosage and PAV
```
for i in `seq 1 30` ; do
index=1
grep -v '^$' chr${i}.group |while IFS= read -r OR; do
printf -v index_padded "%03d" $index
ORid=$(echo R${i}_${index_padded})

####### write coding OR #######
echo $OR |tr ',' '\n'|grep -v '^$' | grep -f - all_coding_annotation.rename.bed.clue |awk -v a=${ORid} '$1=a {print $0}' >> coding.anno

####### write validated pseudo OR #######
echo $OR |tr ',' '\n'|grep -v '^$' | \
grep -f - /home/wgzhang/cactus/bos/new_work/graph_projection/chr${i}.injref.untangle.tsv.collapse | \
awk 'BEGIN{OFS="\t"} {print $1,$3,$4}' | bedtools sort -i - |\
bedtools subtract -A -a - -b /home/wgzhang/cactus/bos/new_work/graph_projection/chr${i}.coding_annotation.bed | \
bedtools merge -i - |\
bedtools intersect -wa -b -  -a /home/wgzhang/cactus/bos/new_work/graph_projection/chr${i}.dead_annotation.bed > validated_pseudo_proto.bed
awk -F'#' -v a=$i '{print $0,$1"_R"a}' validated_pseudo_proto.bed | awk 'BEGIN{OFS="\t"} {print $5"_"$2}' > temp 
cat replace_genome.clue |while IFS= read -r line; do
old=$(echo $line |awk '{print $1}')
new=$(echo $line |awk '{print $2}')
sed -i "s/${old}/${new}/g"  temp
done
awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' validated_pseudo_proto.bed |paste - temp >> validated_pseudo.bed
awk 'BEGIN{OFS="\t"} {print $4}' validated_pseudo_proto.bed |paste - temp >> validated_pseudo.bed.clue
awk -v a=$ORid 'BEGIN{OFS="\t"} {print a,$1}' temp >> validated_pseudo.bed.anno

####### write unvalidated pseudo OR #######
echo $OR |tr ',' '\n'|grep -v '^$' | \
grep -f - /home/wgzhang/cactus/bos/new_work/graph_projection/chr${i}.injref.untangle.tsv.collapse | \
awk 'BEGIN{OFS="\t"} {print $1,$3,$4}' | bedtools sort -i - |\
bedtools subtract -A -a - -b /home/wgzhang/cactus/bos/new_work/graph_projection/chr${i}.coding_annotation.bed | \
bedtools merge -i - |\
bedtools subtract -A -a -  -b /home/wgzhang/cactus/bos/new_work/graph_projection/chr${i}.dead_annotation.bed  > temp 
awk -F '[#\t]'  -v a=$i -v b=$ORid  'BEGIN{OFS="\t"} {print b,$1"_R"a"_"$4}' temp > temp.anno
cat replace_genome.clue |while IFS= read -r line; do
old=$(echo $line |awk '{print $1}')
new=$(echo $line |awk '{print $2}')
sed -i "s/${old}/${new}/g"  temp.anno
done
cut -f2 temp.anno | paste temp - >> unvalidated_pseudo.bed
cat temp.anno >> unvalidated_pseudo.bed.anno

index=$((index + 1))
done
done
```
### 3. heatmap
```
awk '{print $1}' coding.anno > coding.anno.col1
awk '{print $2}' coding.anno | awk -F '_' '{print $1}' |paste coding.anno.col1 - |sort -k1,1 -k2,2 |uniq -c  > coding.count
cat coding.anno unvalidated_pseudo.bed.anno validated_pseudo.bed.anno |awk '{print $1}'  > coding_pseudo.col1
cat coding.anno unvalidated_pseudo.bed.anno validated_pseudo.bed.anno |awk '{print $2}' | awk -F '_' '{print $1}' |paste coding_pseudo.col1 - |sort -k1,1 -k2,2 |uniq -c  > coding_pseudo.count

#R
df=read.table("coding.count",head=F)
or=unique(df$V2)
gen=unique(df$V3)
or_count=as.data.frame(matrix(0,nrow=length(or),ncol=length(gen)))
colnames(or_count)=gen
rownames(or_count)=or
for(i in 1:nrow(df)){
or_count[which(rownames(or_count)==df[i,2]),which(colnames(or_count)==df[i,3])]=df[i,1]
}
or_count <- or_count[apply(or_count, 1, function(x) var(x, na.rm = TRUE)) != 0, ]
sampleanno=read.table("annotation.txt",head=F,sep="\t")
colnames(sampleanno)=c("name","lineage","PC1")
rownames(sampleanno)=sampleanno[,1]
sampleanno[,1]=1
png("PGGB_coding_dosage.png",width=2000,height=2000,res=300)
pheatmap(or_count,annotation_col=sampleanno,color = colorRampPalette(colors = c("white","red"))(100))
dev.off()

setwd("/share/org/YZWL/yzwl_zhangwg/OR_project/1.pangenome/1.3bos_PGGB/pan_cattle31_new/new_work/6.heatmap")
df=read.table("coding_pseudo.count",head=F)
or=unique(df$V2)
gen=unique(df$V3)
or_count=as.data.frame(matrix(0,nrow=length(or),ncol=length(gen)))
colnames(or_count)=gen
rownames(or_count)=or
for(i in 1:nrow(df)){
or_count[which(rownames(or_count)==df[i,2]),which(colnames(or_count)==df[i,3])]=df[i,1]
}
or_count <- or_count[apply(or_count, 1, function(x) var(x, na.rm = TRUE)) != 0, ]
sampleanno=read.table("annotation.txt",head=F,sep="\t")
colnames(sampleanno)=c("name","lineage","PC1")
rownames(sampleanno)=sampleanno[,1]
sampleanno[,1]=1
png("PGGB_coding_pseudo_dosage.png",width=2000,height=2000,res=300)
pheatmap(or_count,annotation_col=sampleanno,color = colorRampPalette(colors = c("white","red"))(100))
dev.off()
```
