nchr=$1


og=chr${nchr}.final.og
$(awk -F' ' -v a=$og '$4==a' ../og_link)
ln -s /share/home/yzwl_zhangwg/OR_project/PGGB/script/coding_annotation_generator.sh ./
ln -s /share/home/yzwl_zhangwg/OR_project/PGGB/script/collaspe_untangle.sh ./
ln -s /share/home/yzwl_zhangwg/OR_project/PGGB/script/dead_annotation_generator.sh ./
ln -s /share/home/yzwl_zhangwg/OR_project/PGGB/script/modify_coor.sh ./
ln -s /share/home/yzwl_zhangwg/OR_project/PGGB/script/OR_counter.sh ./
ln -s /share/home/yzwl_zhangwg/OR_project/PGGB/script/precess_bed.sh ./
ln -s /share/home/yzwl_zhangwg/OR_project/PGGB/script/projection.sh ./
ln -s /share/home/yzwl_zhangwg/OR_project/PGGB/script/pseudo_matrix.sh ./
chmod +x modify_coor.sh
chmod +x precess_bed.sh
chmod +x coding_annotation_generator.sh
chmod +x collaspe_untangle.sh
chmod +x OR_counter.sh
chmod +x dead_annotation_generator.sh
chmod +x projection.sh
chmod +x pseudo_matrix.sh
CHROM=chr${nchr}
fastapath=/share/home/yzwl_zhangwg/OR_project/PGGB/pan_by_chr_new
OR_path=/share/home/yzwl_zhangwg/OR_project/genome2OR
./coding_annotation_generator.sh ${fastapath} ${OR_path} ${CHROM}
./dead_annotation_generator.sh ${fastapath} ${OR_path} ${CHROM}

