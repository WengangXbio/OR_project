source activate odgi
for NC in `cut -f1 ref.info |awk -F"#" '{print $3}'` ; do
grep ${NC} ref.info > ${NC}.bed
odgi extract -i cattle31.full.og -o cattle31.${NC}.full.og -b ${NC}.bed -t 16 -c 0 -E
done

