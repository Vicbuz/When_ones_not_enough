#Bioinformtics steps for WGS processing.
#some scripts were run in loops and some in parallel. Parallel works here are given as an example with one sample for clarity.
#loop work is show as a loop working on a directory of all samples. 

#instead of supplying every single loop and parallel code (as they would just be repeats with different file locations or sample names)
# here a work flow is presented that all samples went through.

#Trimmomatic - trim reads on quality:

for f in /WGS_work/WGS_/*1.fq.gz 
do a=$f
echo $f
IFS=. components=($a)
java -jar trimmomatic-0.39.jar PE -phred33 ${components[0]}.1.fq.gz ${components[0]}.2.fq.gz WGS_work/WGS_trimmo_paired/$(basename $components).1.fq.gz WGS_work/WGS_trimmo_unpaired/$(basename $components).1.fq.gz /WGS_work/WGS_trimmo_paired/$(basename $components).2.fq.gz /WGS_trimmo_unpaired/$(basename $components).2.fq.gz  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
done

#BWA Align to reference genome:

./bwa mem GCF_003254395.2_Amel_HAv3.1_genomic_refseq.fna -M -c 1 /trimmed_paired/116.1.fq.gz /trimmed_paired/116.2.fq.gz > /align_trimmed_paired/116.sam

#Samtools q20 loop:

for f in /WGS_work/WGS_trimmo_paired_align_K/*
do
	echo "Processing $f"
	samtools view -q 20 -h $f > /WGS_work/WGS_trimmo_paired_aligned_q20_K/$(basename "$f")${filename##*.} 
done


#Samtools sort and conversion to bam loop:

for f in /WGS_work/WGS_trimmo_paired_aligned_q20_K/*
do
	echo "Processing $f"
	samtools sort -O bam -o /WGS_work/WGS_trimmo_paired_aligned_q20_sort_K/$(basename "$f")${filename##*.} $f 
done;

#Samtools index:

samtools index -b /WGS_work/WGS_trimmo_paired_aligned_q20_sort_RG_K/116_RG.bam \

#Picard read group per sample:

./gatk AddOrReplaceReadGroups \
I=WGS_work/WGS_trimmo_paired_aligned_q20_sort_K/116.bam \
O=WGS_work/WGS_trimmo_paired_aligned_q20_sort_RG_K/116_RG.bam \
RGID=4 \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=116 

#GATK read in bam to g.vcf per sample:

./gatk HaplotypeCaller \
-I /WGS_work/WGS_trimmo_paired_aligned_q20_sort_RG_K/116_RG.bam \
-R GCF_003254395.2_Amel_HAv3.1_genomic_refseq.fna \
-O /WGS_work/gvcfs/gvcfs_K/116.g.vcf \
-ERC GVCF \

#GATK to per chromosome intermediate database that includes all samples.
# example here is chromosome number 8:

./gatk GenomicsDBImport \
--genomicsdb-workspace-path WGS_work/gvcfs/genomicsdb_eight \
-L WGS_work/gvcfs/eight.list \
--sample-name-map WGS_work/gvcfs/gvcfWGS.sample_map \
--tmp-dir WGS_work/gvcfs/tmp 

#GATK GenotypeGVCFs to make a vcf per chormosome - example here is chromosome number 8.

./gatk GenotypeGVCFs \
-R GCF_003254395.2_Amel_HAv3.1_genomic_refseq.fna \
-V gendb:///WGS_work/gvcfs/genomicsdb_eight \
-O WGS_work/WGSchr8.vcf.gz 

#GATK Gathergvcfs - gather all chromosomes into one vcf.

./gatk GatherVcfs \
I= WGS_work/WGSchr1.vcf.gz \
I= WGS_work/WGSchr2.vcf.gz \
I= WGS_work/WGSchr3.vcf.gz \
I= WGS_work/WGSchr4.vcf.gz \
I= WGS_work/WGSchr5.vcf.gz \
I= WGS_work/WGSchr6.vcf.gz \
I= WGS_work/WGSchr7.vcf.gz \
I= WGS_work/WGSchr8.vcf.gz \
I= WGS_work/WGSchr9.vcf.gz \
I= WGS_work/WGSchr10.vcf.gz \
I= WGS_work/WGSchr11.vcf.gz \
I= WGS_work/WGSchr12.vcf.gz \
I= WGS_work/WGSchr13.vcf.gz \
I= WGS_work/WGSchr14.vcf.gz \
I= WGS_work/WGSchr15.vcf.gz \
I= WGS_work/WGSchr16.vcf.gz \
O= WGS_work/WGS_RAW_nooutgroup.vcf

#GATK SelectVariants - remove non variant sites and select SNPs.

cd programmes/gatk-4.2.0.0/ 
./gatk SelectVariants \
-R GCF_003254395.2_Amel_HAv3.1_genomic_refseq.fna \
-V WGS_work/WGS_RAW_nooutgroup.vcf \
--select-type-to-include SNP \
--exclude-non-variants true \
-O WGS_work/WGS_RAW_nooutgroup_SNPs.vcf

#GATK VariantFiltration - Filter using the GATK recommended thresholds.

./gatk VariantFiltration \
-V WGS_work/WGS_RAW_nooutgroup_SNPs.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered.vcf


#vctools - filter out sites with exsessive depth. Get report from --site-depth and visualise in R or python. 

vcftools --vcf WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered.vcf --site-depth --out WGS_RAW_nooutgroup_SNPs_gatkfiltered

#vcftools filter out min and max depth

vcftools --vcf WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered.vcf --min-alleles 2 \
--max-alleles 2 \
--max-meanDP 65 \
--min-meanDP 10 \
--recode --recode-INFO-all --out WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2

#vcftools - genotypes need to be in 90% of samples

vcftools --vcf WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2.recode.vcf \
--max-missing 0.9 \
--recode --recode-INFO-all --out WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90 

# makes a missingness file

vcftools --vcf WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_maxmissing_90.vcf.recode.vcf \
--missing-indv --out WGS_work/missingness_WGS 

# cutstom python code to find sample with % genotype - individuals needs to have 80% of all genotypes

python imiss_rmv_samples.py WGS_work/missingness_WGS.imiss 0.8 WGS_work/WGSto_be_rmv_80.txt 

#individuals needs to have 80% of all genotypes

vcftools --vcf WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90.recode.vcf \
--remove WGS_work/WGSto_be_rmv_80.txt \
--recode --recode-INFO-all --out WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90_indvidmissing80


