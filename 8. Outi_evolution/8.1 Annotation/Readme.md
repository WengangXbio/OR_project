### 1. run genome2or
```
cd /share/org/YZWL/yzwl_zhangwg/OR_project/2.outi_evolution/1.annotation/g2or/
sh run_g2or.sh
```

### 2. annotations
```
cd /share/org/YZWL/yzwl_zhangwg/OR_project/2.outi_evolution/1.annotation/annotation/

g2or="/share/org/YZWL/yzwl_zhangwg/OR_project/2.outi_evolution/1.annotation/g2or"
for genome in `cat /share/org/YZWL/yzwl_zhangwg/OR_project/2.outi_evolution/1.annotation/annotation/select_genome` ; do
cat ${g2or}/${genome}.g2or/ORannotation_itera1_final_low-quality.fasta ${g2or}/${genome}.g2or/ORannotation_itera1_final_pseu_ORs.fasta > ${genome}.dead.fasta        
cp ${g2or}/${genome}.g2or/ORannotation_itera2_final_func_dna_ORs.fasta ${genome}.func_dna_ORs.fasta
cp ${g2or}/${genome}.g2or/ORannotation_itera2_final_func_pro_ORs.fasta  ${genome}.func_pro_ORs.fasta
done

cat ../select_genome_ABR | while IFS= read -r line; do
genome=$(echo $line |awk '{print $1}')
name=$(echo  $line |awk '{print $2}')
samtools faidx ${genome}.func_pro_ORs.fasta
awk -v a=${name} '{printf a"-""%04d %s\n",NR,$1}' ${genome}.func_pro_ORs.fasta.fai |awk 'BEGIN{OFS="\t"} {print $2,$1}' > ${genome}.func_pro_ORs.fasta.rename.trace
seqkit replace -p "(\S+)" -r "{kv}" -k${genome}.func_pro_ORs.fasta.rename.trace ${genome}.func_pro_ORs.fasta > ${genome}.func_pro_ORs.fasta.rename.fasta
samtools faidx ${genome}.func_dna_ORs.fasta
seqkit replace -p "(\S+)" -r "{kv}" -k${genome}.func_pro_ORs.fasta.rename.trace ${genome}.func_dna_ORs.fasta > ${genome}.func_dna_ORs.fasta.rename.fasta
samtools faidx ${genome}.dead.fasta
awk -v a=${name} '{printf a"-dead-""%04d %s\n",NR,$1}' ${genome}.dead.fasta.fai |awk 'BEGIN{OFS="\t"} {print $2,$1}' > ${genome}.dead.fasta.rename.trace
seqkit replace -p "(\S+)" -r "{kv}" -k${genome}.dead.fasta.rename.trace ${genome}.dead.fasta > ${genome}.dead.fasta.rename.fasta
done

cat *.func_pro_ORs.fasta.rename.fasta > all_36outi.func_pro_ORs.fasta.rename.fasta
cat *.func_dna_ORs.fasta.rename.fasta > all_36outi.func_dna_ORs.fasta.rename.fasta
cat *.dead.fasta.rename.fasta > all_36outi.dead.fasta.rename.fasta
```

