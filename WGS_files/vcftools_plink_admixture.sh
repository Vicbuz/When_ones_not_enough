#!/bin/bash

#SBATCH -J vcffilter               # Job name
#SBATCH --partition normal               # Job queue
#SBATCH -o job.%j.out         # Name of stdout output file (%j expands to jobId)
#SBATCH -N 1                  # Total number of nodes requested
#SBATCH -n 32                 # Total number of mpi tasks requested
#SBATCH -t 72:00:00           # Run time (hh:mm:ss) - 1.5 hours


#code will remove the buckfast samples and replicated that are not needed in final analysis

vcftools --vcf /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90_indvidmissing80.recode.vcf \
--remove /home/vbuswell/WGS_work/rep_buck.txt \
--recode --recode-INFO-all --out /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90_indvidmissing80_noreps_nobuck

#arbitary thin the data 1kb apart

vcftools --vcf /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90_indvidmissing80_noreps_nobuck.recode.vcf \
--thin 1000 \
--recode --recode-INFO-all --out /home/vbuswell/WGS_work/wgs_1kb_british_isles_only.recode.vcf

#plink doesn't accept chromosome numbers must be renamed i.e. chr1 

bcftools annotate --rename-chrs /home/vbuswell/WGS_work/chr_rename_map.txt /home/vbuswell/WGS_work/wgs_1kb_british_isles_only.recode.vcf \
-Ov -o /home/vbuswell/WGS_work/chr_edited_WGS_1kb_thinned.vcf

#make a set of plink files for admixture:

vcftools --vcf /home/vbuswell/WGS_work/chr_edited_WGS_1kb_thinned.vcf --plink --out /home/vbuswell/WGS_work/WGS_1kb_thinned

export PATH=$PATH/home/vbuswell/programmes/

cd /home/vbuswell/programmes/

#make a bed file for admixture

./plink --file /home/vbuswell/WGS_work/WGS_1kb_thinned  --make-bed --out /home/vbuswell/WGS_work/plink_final_filtered_WGS_1kb_thinned --chr-set 32 no-xy --allow-extra-chr

#run adimxture

cd /home/vbuswell/admixture_linux-1.3.0

for K in 1 2 3 4 5 6 7 8 9 10; \
do ./admixture --cv /home/vbuswell/WGS_work/plink_final_filtered_WGS_1kb_thinned.bed $K | tee log${K}.out;  #this will crossvalidate K to indicate the most appropriate K for the data
done

grep -h CV log*.out