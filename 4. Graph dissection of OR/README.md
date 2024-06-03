## Graph dissection of OR genes
### PGGB-based graph
#### To dissect graph of OR regions, there are 4 steps
#### 
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



