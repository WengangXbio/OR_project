awk '{print $3}' ../combined_ORgt.dosage |sort |uniq > combined_ORgt.OR.name
for id in `ls *vcf.or |awk -F'.' '{print $1}'`; do
echo -e "or ${id}" |cat - combined_ORgt.OR.name |awk '{print $0, 0}' > ${id}.dosage
cat ${id}.PanGenie_genotyping.vcf.or |while IFS= read -r line; do
orid=$(echo $line |awk '{print $1"_"$2}')
gt1=$(echo $line |awk '{print $3}' |awk -F'/' '{print "GT"$1}')
gt2=$(echo $line |awk '{print $3}' |awk -F'/' '{print "GT"$2}')
if [ "$gt1" == "GT." ]; then
awk -v a=$orid '{if($1 ~ "^"a) $2="NA"; print}' ${id}.dosage > temp0 && mv temp0 ${id}.dosage
else
awk -v a=$orid -v b=$gt1 '$2==a && $1==b' ../combined_ORgt.dosage |while IFS= read -r vg1; do
tid=$(echo $vg1 |awk '{print $3}')
dog=$(echo $vg1 |awk '{print $4}')
awk -v a="$tid" -v b="$dog" '{if ($1==a) print $1, $2+b; else print $0}' ${id}.dosage > temp
mv temp ${id}.dosage
done
awk -v a=$orid -v b=$gt2 '$2==a && $1==b' ../combined_ORgt.dosage|while IFS= read -r vg1; do
tid=$(echo $vg1 |awk '{print $3}')
dog=$(echo $vg1 |awk '{print $4}')
awk -v a="$tid" -v b="$dog" '{if ($1==a) print $1, $2+b; else print $0}' ${id}.dosage > temp
mv temp ${id}.dosage
done
fi
done
done


for file in *.dosage; do
awk '{print $2}' "$file" > "${file}.tmp"
done
paste *.tmp  > combined.txt
rm *.tmp
