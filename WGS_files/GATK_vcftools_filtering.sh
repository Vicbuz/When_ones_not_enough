#!/bin/bash

#SBATCH -J Gth.Slct.flt               # Job name
#SBATCH --partition normal               # Job queue
#SBATCH -o job.%j.out         # Name of stdout output file (%j expands to jobId)
#SBATCH -N 1                  # Total number of nodes requested
#SBATCH -n 32                 # Total number of mpi tasks requested
#SBATCH -t 72:00:00           # Run time (hh:mm:ss) - 1.5 hours

# Launch MPI-based executable
export LD_LIBRARY_PATH=/home/vbuswell/anaconda3/lib:$LD_LIBRARY_PATH

export PATH=$PATH/home/vbuswell/programmes/gatk-4.1.9.0
export PATH=$PATH:/home/vbuswell/anaconda3/pkgs/vcftools-0.1.14-pl5.22.0_0/bin

cd /home/vbuswell/programmes/gatk-4.2.0.0/ 
./gatk GatherVcfs \
I= /home/vbuswell/WGS_work/WGSchr1.vcf.gz \
I= /home/vbuswell/WGS_work/WGSchr2.vcf.gz \
I= /home/vbuswell/WGS_work/WGSchr3.vcf.gz \
I= /home/vbuswell/WGS_work/WGSchr4.vcf.gz \
I= /home/vbuswell/WGS_work/WGSchr5.vcf.gz \
I= /home/vbuswell/WGS_work/WGSchr6.vcf.gz \
I= /home/vbuswell/WGS_work/WGSchr7.vcf.gz \
I= /home/vbuswell/WGS_work/WGSchr8.vcf.gz \
I= /home/vbuswell/WGS_work/WGSchr9.vcf.gz \
I= /home/vbuswell/WGS_work/WGSchr10.vcf.gz \
I= /home/vbuswell/WGS_work/WGSchr11.vcf.gz \
I= /home/vbuswell/WGS_work/WGSchr12.vcf.gz \
I= /home/vbuswell/WGS_work/WGSchr13.vcf.gz \
I= /home/vbuswell/WGS_work/WGSchr14.vcf.gz \
I= /home/vbuswell/WGS_work/WGSchr15.vcf.gz \
I= /home/vbuswell/WGS_work/WGSchr16.vcf.gz \
O= /home/vbuswell/WGS_work/WGS_RAW_nooutgroup.vcf

cd /home/vbuswell/programmes/gatk-4.2.0.0/ 
./gatk SelectVariants \
-R /home/vbuswell/novaseq/gatk/GCF_003254395.2_Amel_HAv3.1_genomic_refseq.fna \
-V /home/vbuswell/WGS_work/WGS_RAW_nooutgroup.vcf \
--select-type-to-include SNP \
--exclude-non-variants true \
-O /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs.vcf

#Here we show basic filtering thresholds researchers may find useful to start.

cd /home/vbuswell/programmes/gatk-4.2.0.0/ 

./gatk VariantFiltration \
-V /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered.vcf

# cd /home/vbuswell/WGS_work/

vcftools --vcf /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered.vcf --min-alleles 2 --max-alleles 2 --max-meanDP 65 --min-meanDP 10 --recode --recode-INFO-all --out /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2

vcftools --vcf /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2.recode.vcf --max-missing 0.9 --recode --recode-INFO-all --out /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90 #genotypes need to be in 90% of samples

vcftools --vcf /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_maxmissing_90.vcf.recode.vcf --missing-indv --out /home/vbuswell/WGS_work/missingness_WGS # makes a missingness file

python imiss_rmv_samples.py /home/vbuswell/WGS_work/missingness_WGS.imiss 0.8 /home/vbuswell/WGS_work/WGSto_be_rmv_80.txt # individuals needs to have 80% of all genotypes

vcftools --vcf /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90.recode.vcf --remove /home/vbuswell/WGS_work/WGSto_be_rmv_80.txt --recode --recode-INFO-all --out /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90_indvidmissing80

cd /home/vbuswell/WGS_work/

vcftools --vcf /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90_indvidmissing80 
