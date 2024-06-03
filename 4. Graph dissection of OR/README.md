## Graph dissection of OR genes
### PGGB-based graph
#### To dissect graph of OR regions, there are 4 steps
#### Processing OR annotations by adjusting coordinates, and create projection of OR among assemblies
```
for nchr in `seq 1 29` ; do
mkdir chr${nchr} -p
cd chr${nchr}
og=chr${nchr}.final.og
ln -s /share/home/yzwl_zhangwg/OR_project/PGGB/script/preperation.sh ./
./preperation.sh ${nchr}
sed -i "s/OMG/${og}/g" projection.sh
csub < projection.sh
cd ..
done
```
#### Adjusting OR projections by collapsing coordinates
```
for nchr in `seq 1 29` ; do
cd chr${nchr}
if [ -s "./injref.untangle.tsv" ]; then
cp /share/home/yzwl_zhangwg/OR_project/PGGB/script/count_run.sh ./
sed -i "s/CHRN/${nchr}/g" count_run.sh
chmod +x count_run.sh
csub < count_run.sh
else 
   echo "chr${nchr} is blank"
fi
cd ..
done
```
#### Counting OR dosages and living or dead
```
for nchr in `seq 1 29` ; do
echo chr${nchr}
cd chr${nchr}
if [ -s "./injref.untangle.tsv.collapse" ]; then
cp -f /share/home/yzwl_zhangwg/OR_project/PGGB/script/count.matrix.sh ./
cp -f /share/home/yzwl_zhangwg/OR_project/PGGB/script/pseudo_matrix.sh ./
cp -f /share/home/yzwl_zhangwg/OR_project/PGGB/script/OR_counter.sh ./
sh count.matrix.sh ${nchr}
else 
   echo "chr${nchr} is blank"
fi
cd ..
done
```

#### Collection of OR information from assemblies
```
cat genomelist | tr '\n' '\t' > or.counts.matrix0.collection
cat genomelist | tr '\n' '\t' > or.pseudo.matrix0.collection
printf '\n' >> or.counts.matrix0.collection
printf '\n' >> or.pseudo.matrix0.collection
for nchr in `seq 1 29` ; do
cd chr${nchr}
if [ -s "./or.counts.matrix0" ]; then
sed '1d' or.counts.matrix0 |cat ../or.counts.matrix0.collection - > ../temp
mv ../temp ../or.counts.matrix0.collection
sed '1d' or.pseudo.matrix0 |cat ../or.pseudo.matrix0.collection - > ../temp
mv ../temp ../or.pseudo.matrix0.collection
else 
   echo "chr${nchr} is blank"
fi
cd ..
done
```
#### Collapsing coordinates is slow when OR is highly enriched in certain chromosome (e.g. chromosome 15), running large collasping scripts as follow
```
for iseq in `cat inputseq` ; do
mkdir ${iseq}_collaspe 
cd ${iseq}_collaspe 
cp /share/home/yzwl_zhangwg/OR_project/PGGB/script/collaspe_untangle.sh collaspe_untangle.sh 
cp /share/home/yzwl_zhangwg/OR_project/PGGB/script/large_collaspe.sh large_collaspe.sh
cp ../injref.untangle.tsv ./
echo $iseq > inputseq
csub < large_collaspe.sh
cd ..
done
```
```
echo > collection 
for splitfile in `ls -l ./ | grep '^d' | grep "_collaspe" | awk '{print $9}'`; do
cat collection ${splitfile}/injref.untangle.tsv.collapse > temp
mv temp collection
done
mv collection injref.untangle.tsv.collapse
```




