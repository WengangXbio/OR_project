## OR annotation by genome2OR
```
cd /share/home/yzwl_zhangwg/OR_project/genome2OR
for genome in `ls *.fna`  ; do
outdir="${genome}.g2or"
sed "s/INPUT/${genome}/g" genome2OR.sh | sed "s/OUTPUT/${outdir}/g" > ${genome}.sh
csub < ${genome}.sh
done
```
```
[yzwl_zhangwg@mgt15 genome2OR]$ cat genome2OR.sh 
#!/bin/bash
#CSUB -J genome2OR
#CSUB -q cpu
#CSUB -o %J.out
#CSUB -e %J.error
#CSUB -n 24   #core
#CSUB -R span[hosts=1]

module load anaconda3/4.12.0
source activate genome2or
Iteration Mammalia OUTPUT INPUT 
```
