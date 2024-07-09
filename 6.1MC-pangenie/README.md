## MC + Pangenie workflow && dissection bubble with OR genes 
### Note
1. MC graph construction can be found [here](https://github.com/WengangXbio/OR_project/blob/main/6.1MC-pangenie/scripts/MC_construction.sh)
2. MC produced a raw vcf with **nested bubble** and **haploid genotyping**, thus should be modified as overlapping variants with vcfbub and to be diploid with [shell scripts](https://github.com/WengangXbio/OR_project/blob/main/6.1MC-pangenie/scripts/vcf_modify.sh).
3. Removing alternative genotype with sequence containing "N", using the [script here](https://github.com/WengangXbio/OR_project/blob/main/6.1MC-pangenie/scripts/removeN_vcf.sh)
4. Modified vcf is indexed and run yazhouwan sever ([see scripts here](https://github.com/WengangXbio/OR_project/blob/main/6.1MC-pangenie/scripts/pangenie-index.sh)) and local server ([see scripts here](https://github.com/WengangXbio/OR_project/blob/main/6.1MC-pangenie/scripts/pangenie-index_local.sh))
5. Find bubbles with OR genes with [script](https://github.com/WengangXbio/OR_project/blob/main/6.1MC-pangenie/scripts/find_OR_bubble.sh);
6. Dissection OR bubbles including finding ortholog OR genes among alternative path and counting OR genes for each alternative path [see script here](https://github.com/WengangXbio/OR_project/blob/main/6.1MC-pangenie/scripts/OR_dosage_and_genotype.sh)
7. Extract OR bubbles from pangenie output vcf ([script](https://github.com/WengangXbio/OR_project/blob/main/6.1MC-pangenie/scripts/OR_bubble_vcf_extraction.sh))
8. 

