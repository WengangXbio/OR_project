## PGGB pipeline

PGGB is implemented by one-by-one chromosome, the first step is to decide which chromosomes/scaffolds/cotigs are corresponding to targeted chromosome (ARS2.0) for each assemblies.

Here, a script called KMAP is written for sorting chromosomes/scaffolds/cotigs.

### Step1. Build consolidate alignment map corresponding to reference for each assembly
### Step2. Extract chromosomes/scaffolds/cotigs from assemblies for chromosome-pangenome
### Step3. Submit PGGB job in batch

<img src="https://github.com/WengangXbio/OR_project/blob/f4667915ec81e06fa33b267a95dfe55b8c876e0d/3.%20PGGB/Figures/kmap_principal.png">

#### Run step 1, built alignment map
```
cd /share/home/yzwl_zhangwg/OR_project/PGGB/assembly/NGenomeSyn
for fna in `ls *.paf |awk -F'.paf' '{print $1}'` ; do
sed "s/ASSEMBLY/$fna/g" kmap_run.sh > kmap_run$fna.sh
csub < kmap_run$fna.sh
done
```
#### Run step 2, extract coordinates from assemblies by chromosome
```
for i in `seq 1 31`  ; do
chrn=chr$i
sed "s/CHRN/$chrn/g" kmap_run2.sh > kmap_run2_$chrn.sh
csub < kmap_run2_$chrn.sh
done
```
