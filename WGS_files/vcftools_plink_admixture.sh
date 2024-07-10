#!/bin/bash

#code will remove the buckfast samples and replicated that are not needed in final analysis

vcftools --vcf WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90_indvidmissing80.recode.vcf \
--remove WGS_work/rep_buck.txt \
--recode --recode-INFO-all --out WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90_indvidmissing80_noreps_nobuck

#arbitary thin the data 1kb apart

vcftools --vcf WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90_indvidmissing80_noreps_nobuck.recode.vcf \
--thin 1000 \
--recode --recode-INFO-all --out WGS_work/wgs_1kb_british_isles_only.recode.vcf

#plink doesn't accept chromosome numbers must be renamed i.e. chr1 

bcftools annotate --rename-chrs WGS_work/chr_rename_map.txt /home/vbuswell/WGS_work/wgs_1kb_british_isles_only.recode.vcf \
-Ov -o WGS_work/chr_edited_WGS_1kb_thinned.vcf

#make a set of plink files for admixture:

vcftools --vcf WGS_work/chr_edited_WGS_1kb_thinned.vcf --plink --out WGS_work/WGS_1kb_thinned

#make a bed file for admixture

./plink --file WGS_work/WGS_1kb_thinned  --make-bed --out WGS_work/plink_final_filtered_WGS_1kb_thinned --chr-set 32 no-xy --allow-extra-chr

#run adimxture

cd admixture_linux-1.3.0

for K in 1 2 3 4 5 6 7 8 9 10; \
do ./admixture --cv WGS_work/plink_final_filtered_WGS_1kb_thinned.bed $K | tee log${K}.out;  #this will crossvalidate K to indicate the most appropriate K for the data
done

grep -h CV log*.out
