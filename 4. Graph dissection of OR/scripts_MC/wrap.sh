NC=NNNCCC
nchr=NCHRNUM
cp -f ../collaspe_untangle.sh ./
cp -f ../cattle31.${NC}.full.og.L  inputseq
chmod +x collaspe_untangle.sh
./collaspe_untangle.sh injref.untangle.tsv inputseq injref.untangle.tsv.collapse
cp -f ../count.matrix.sh ./
cp -f ../pseudo_matrix.sh ./
cp -f ../OR_counter.sh ./
chmod +x count.matrix.sh 
chmod +x pseudo_matrix.sh
chmod +x OR_counter.sh
sh count.matrix.sh ${nchr} ${NC}

