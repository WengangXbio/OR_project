cat or_anno_transfer/*.func.bed > all.func.bed
cat or_anno_transfer/*.dead.bed > all.dead.bed
awk -F'#0' '{ if ($1 == "GCF_002263795") {print $1"#0"$2, $3} else {print $0} }' all.func.bed |awk -F' ' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' > all.func.bed0
awk -F'#0' '{ if ($1 == "GCF_002263795") {print $1"#0"$2, $3} else {print $0} }' all.dead.bed |awk -F' ' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' > all.dead.bed0

for NC in `cut -f1 ref.info |awk -F"#" '{print $3}'` ; do
awk -F':' '{print $1}' cattle31.${NC}.full.og.L |grep -f - all.func.bed0 > ${NC}.func.bed
awk -F':' '{print $1}' cattle31.${NC}.full.og.L |grep -f - all.dead.bed0 > ${NC}.dead.bed
odgi procbed -i cattle31.${NC}.full.og -b ${NC}.func.bed -t 4 > ${NC}.func.proc.bed
odgi procbed -i cattle31.${NC}.full.og -b ${NC}.dead.bed -t 4 > ${NC}.dead.proc.bed
done

