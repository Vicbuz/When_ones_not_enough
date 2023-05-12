# When_ones_not_enough
Code used in paper:  When one’s not enough: Colony pool-seq outperforms individual-based methods for assessing introgression in Apis mellifera subspecies

https://www.mdpi.com/2075-4450/14/5/421 

csv_to_ped.py was used to convert the SNP array raw data into a ped file for use with ADMIXTURE.

vcf_to_colony_allele_frequencies.py takes an input vcf file created by GATK, calculates allele frequencies from the AD column and outputs a large dataframe of allele frequencies for all samples. This scrip is envoked by:
```bash
python vcf_to_colony_allele_frequencies.py file_in, file_out, snp_num
```

abba_baba_colony.py perfroms D statistics, significance testing and f statistics on colony allele frequencies data frame obtained from a vcf file format. The input data frame for which needs to be sample names as a header and rows of SNPs.
