for i in `seq 1 31`  ; do
chr="chr${i}"
NC=$(grep $chr ARS20_chr_NC.info|awk '{print $2}')
cd $chr
sed "s/CHRNUM/${chr}/g" ../pggb.sh  > pggb.$chr.sh
csub < pggb.$chr.sh
cd ..
done

