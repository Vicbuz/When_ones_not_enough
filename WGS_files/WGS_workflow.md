# Bioinformtics workflow for WGS processing.

This workflow is the detailed methods for the paper Buswell et al 2024, Whole genome analyses of introgression in British and Irish Apis mellifera mellifera, in review DOI incoming.

Some scripts were run in loops and some in parallel. Parallel works here are given as an example with one sample for clarity, loops remain.

This workflow includes the following software:

Trimmomatic v. 0.39\
BWA MEM v. 0.7.17\
Samtools v. 1.10\
GATK 4.1.9.0\
vcftools v. 0.1.14\
Plink v. 1.07\
ADMIXTURE v. 1.3.0

Notes: follows the GATK4 Best Practice Workflow, using GATK version 4.1.9.0 and was written in 2021. GATK advice is that recommendations evolve in step with the rapid pace of technological and methodological innovation in the field of genomics. i.e. best practices change with time, please always check the up to date best practices. 

If you use any of the above tools you must cite the papers and version of those tools. The authors of this repository are not responsible for the above tools.  

If you use any custom python code from this repository, in the spirit of open science, you may do so freely. We simply request you cite the paper it was originally used in either: Buswell et al 2023 [doi](https://doi.org/10.3390/insects14050421) or Buswell et al 2024 doi incoming (in review).

## Trimmomatic

Files were trimmed using Trimmomatic a loop for every forward and reverse read in a directory.
```
for f in /WGS_work/WGS_/*1.fq.gz 
do a=$f
echo $f
IFS=. components=($a)
java -jar trimmomatic-0.39.jar PE -phred33 ${components[0]}.1.fq.gz ${components[0]}.2.fq.gz WGS_work/WGS_trimmo_paired/$(basename $components).1.fq.gz WGS_work/WGS_trimmo_unpaired/$(basename $components).1.fq.gz /WGS_work/WGS_trimmo_paired/$(basename $components).2.fq.gz /WGS_trimmo_unpaired/$(basename $components).2.fq.gz  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
done
```
## BWA MEM align to reference genome:

Samples were aligned to the genome using bwa in parallel, the example here is for each sample. See [bwa manual for more details](https://bio-bwa.sourceforge.net/bwa.shtml)
```
./bwa mem GCF_003254395.2_Amel_HAv3.1_genomic_refseq.fna -M -c 1 /trimmed_paired/116.1.fq.gz /trimmed_paired/116.2.fq.gz > /align_trimmed_paired/116.sam
```
## Samtools q20 :

Samtools was used in a loop to filter based on mapping quality. See [Samtools manual for more details](http://www.htslib.org/doc/samtools.html)
```
for f in /WGS_work/WGS_trimmo_paired_align_K/*
do
	echo "Processing $f"
	samtools view -q 20 -h $f > /WGS_work/WGS_trimmo_paired_aligned_q20_K/$(basename "$f")${filename##*.} 
done
```

## Samtools sort

Samtools was used to sort by genomic coordinates and output as a bam file. See [Samtools manual for more details](http://www.htslib.org/doc/samtools.html)
```
for f in /WGS_work/WGS_trimmo_paired_aligned_q20_K/*
do
	echo "Processing $f"
	samtools sort -O bam -o /WGS_work/WGS_trimmo_paired_aligned_q20_sort_K/$(basename "$f")${filename##*.} $f 
done;
```
## Samtools index

Samtools was used to index the genomic coordinates for fast random access. See [Samtools manual for more details](http://www.htslib.org/doc/samtools.html)

```
samtools index -b /WGS_work/WGS_trimmo_paired_aligned_q20_sort_RG_K/116_RG.bam \
```
## Picard read group

PICARD tools (available as part of the GATK package) was used to edit read groups per sample. For options and guidance see [GATK AddOrReplaceReadGroups (Picard)](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard)
```
./gatk AddOrReplaceReadGroups \
I=WGS_work/WGS_trimmo_paired_aligned_q20_sort_K/116.bam \
O=WGS_work/WGS_trimmo_paired_aligned_q20_sort_RG_K/116_RG.bam \
RGID=4 \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=116 
```
## GATK read in bam to g.vcf

Bam files were then read into GATK to create a per sample g.vcf. For more details and options see [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360050814612-HaplotypeCaller)
```
./gatk HaplotypeCaller \
-I /WGS_work/WGS_trimmo_paired_aligned_q20_sort_RG_K/116_RG.bam \
-R GCF_003254395.2_Amel_HAv3.1_genomic_refseq.fna \
-O /WGS_work/gvcfs/gvcfs_K/116.g.vcf \
-ERC GVCF \
```
## GATK intermediate database 
This was a run to import single-sample GVCFs into GenomicsDB before joint genotyping. This was run per chromosome and the example here is chromosome number 8. For more details and options see [GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360051305591-GenomicsDBImport)
```
./gatk GenomicsDBImport \
--genomicsdb-workspace-path WGS_work/gvcfs/genomicsdb_eight \
-L WGS_work/gvcfs/eight.list \
--sample-name-map WGS_work/gvcfs/gvcfWGS.sample_map \
--tmp-dir WGS_work/gvcfs/tmp 
```
## GATK GenotypeGVCFs 

From the database a vcf is produced containing all files per chromosome (example here is chromosome number 8). For more detail see [GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360050816072-GenotypeGVCFs)
```
./gatk GenotypeGVCFs \
-R GCF_003254395.2_Amel_HAv3.1_genomic_refseq.fna \
-V gendb:///WGS_work/gvcfs/genomicsdb_eight \
-O WGS_work/WGSchr8.vcf.gz 
```
## GATK Gathergvcfs

This step gathers all chromosomes into one vcf. For more detail see [GatherVcfs (Picard) ](https://gatk.broadinstitute.org/hc/en-us/articles/360050814232-GatherVcfs-Picard)
```
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
```
## GATK SelectVariants 

Select variant sites and SNPs for downstream work. For more guidance see [GATK SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360051305531-SelectVariants)
```
cd programmes/gatk-4.2.0.0/ 
./gatk SelectVariants \
-R GCF_003254395.2_Amel_HAv3.1_genomic_refseq.fna \
-V WGS_work/WGS_RAW_nooutgroup.vcf \
--select-type-to-include SNP \
--exclude-non-variants true \
-O WGS_work/WGS_RAW_nooutgroup_SNPs.vcf
```
# GATK VariantFiltration

Hard-filtering variant calls based on certain criteria using the GATK recommended thresholds. For more detail see [VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360050815032-VariantFiltration)
```
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
```

## vctools depth

Filter out sites with excessive depth indicative of repeat regions. Get report from --site-depth and visualise in R or python to make informed decision. 
```
vcftools --vcf WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered.vcf --site-depth --out WGS_RAW_nooutgroup_SNPs_gatkfiltered
```

Minimum depth cutoffs will remove false positive calls. A maximum depth cut off is important because regions with high read depths are likely repetitive ones with the potential to belong to multiple parts of the genome. 
```
vcftools --vcf WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered.vcf --min-alleles 2 \
--max-alleles 2 \
--max-meanDP 65 \
--min-meanDP 10 \
--recode --recode-INFO-all --out WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2
```
## vcftools missingness

To allow comparison across samples with ADMIXTURE, genotypes need to be in 90% of samples
```
vcftools --vcf WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2.recode.vcf \
--max-missing 0.9 \
--recode --recode-INFO-all --out WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90 
```
vcftools report missingness values of samples
```
vcftools --vcf WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_maxmissing_90.vcf.recode.vcf \
--missing-indv --out WGS_work/missingness_WGS 
```
custom python code to grab names of samples that have less than a given % of genotypes - individuals in the example below needs to have 80% of all genotypes (0.8)
```
python imiss_rmv_samples.py WGS_work/missingness_WGS.imiss 0.8 WGS_work/WGSto_be_rmv_80.txt 
```
Use vcftools and the text file written out by the individuals needs to have 80% of all genotypes.
```
vcftools --vcf WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90.recode.vcf \
--remove WGS_work/WGSto_be_rmv_80.txt \
--recode --recode-INFO-all --out WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90_indvidmissing80
```
## Filter samples and thin data

code will remove the buckfast samples and replicated that are not needed in final analysis
```
vcftools --vcf /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90_indvidmissing80.recode.vcf \
--remove /home/vbuswell/WGS_work/rep_buck.txt \
--recode --recode-INFO-all --out /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90_indvidmissing80_noreps_nobuck
```
This command will arbitary thin the data 1kb apart. 

vcftools --vcf /home/vbuswell/WGS_work/WGS_RAW_nooutgroup_SNPs_gatkfiltered_depth10_65_alleles2_maxmissing90_indvidmissing80_noreps_nobuck.recode.vcf \
--thin 1000 \
--recode --recode-INFO-all --out /home/vbuswell/WGS_work/wgs_1kb_british_isles_only.recode.vcf

## Plink rename files

Plink doesn't accept chromosome naming convention of our data, these were therefore renamed.
```
bcftools annotate --rename-chrs /home/vbuswell/WGS_work/chr_rename_map.txt /home/vbuswell/WGS_work/wgs_1kb_british_isles_only.recode.vcf \
-Ov -o /home/vbuswell/WGS_work/chr_edited_WGS_1kb_thinned.vcf
```
## Make plink files

make a set of plink files for admixture:
```
vcftools --vcf /home/vbuswell/WGS_work/chr_edited_WGS_1kb_thinned.vcf --plink --out /home/vbuswell/WGS_work/WGS_1kb_thinned
```
make a bed file for admixture
```
plink --file /home/vbuswell/WGS_work/WGS_1kb_thinned  --make-bed --out /home/vbuswell/WGS_work/plink_final_filtered_WGS_1kb_thinned --chr-set 32 no-xy --allow-extra-chr
```
## ADMIXTURE 

Run admxiture from K=1 to K=10 and cross validate K to indicate the most appropriate K for the data

```
for K in 1 2 3 4 5 6 7 8 9 10; \
do ./admixture --cv /home/vbuswell/WGS_work/plink_final_filtered_WGS_1kb_thinned.bed $K | tee log${K}.out;  
done

grep -h CV log*.out
```
