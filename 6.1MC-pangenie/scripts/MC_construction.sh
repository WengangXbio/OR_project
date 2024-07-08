source ~/software/cactus-bin-v2.8.2/cactus_env/bin/activate

bsub -J cactus -n 48 -R "span[hosts=1] rusage[mem=10GB] select[maxmem>480G]" -o %J.out -e %J.err -q smp \
"
cactus-pangenome ./js ./genomelist --reference GCF_002263795 --outDir cattle31 --outName cattle31 \
--vcf --giraffe --gfa --gbz --odgi --workDir ./work \
--maxDisk 2000G --maxMemory 480G --maxCores 48 --mgCores 48 --mapCores 48 --consCores 48 \
--consMemory 480G --indexCores 48 --indexMemory 480G
"
