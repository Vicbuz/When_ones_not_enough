# When one's not enough
Code used in paper:  When one’s not enough: Colony pool-seq outperforms individual-based methods for assessing introgression in Apis mellifera subspecies foudn here https://www.mdpi.com/2075-4450/14/5/421 

All code is written in python 3.

### SNP Array format to ped file
The SNP Array data was recieved in a csv format. csv_to_ped.py was used to convert the SNP array raw data into a ped file for use with ADMIXTURE. The use of this function in the script is like so:

```
iPLEX_to_ped('/path/to/your/csv.csv', '/out/path/to/your/new/ped.txt')
```

### Obtain colony allele frequecnies from VCF

vcf_to_colony_allele_frequencies.py takes an input vcf file created by GATK, calculates allele frequencies from the AD and DP fields in the info column and outputs a large dataframe of allele frequencies for all samples. This scrip is envoked by:
```bash
python vcf_to_colony_allele_frequencies.py file_in.vcf file_out.txt snp_num
```
Where ```file_in.vcf``` is your filtered vcf file, ```file_out.txt``` will be the resulting allele freuqencie table and the ```snp_num``` is the number of SNPs in your filtered vcf input. The total number of SNPs is used to as a sense check against the resulting the length of the allele frequency table. The resulting table will be a data frame with sample names as headers and rows of SNPs containing allele frequecnies.

### Perform ABBA BABA

abba_baba_colony.py perfroms D statistics, significance testing and f statistics on colony allele frequencies data frame obtained from a vcf file format (which can be created by using the vcf_to_colony_allele_frequencies.py). The input data frame for which needs to be sample names as a header and rows of SNPs.

The top of the code requires you to place your samples for testing into the lists. Then at the bottom of the code, uncomment the loop procedures you wish to perform. 

WARNING, this code can produce a lot of out files and take a long time to run as each colony is tested multiple times. Work is underway to optimise this code and prevent the long running time and multiple printouts. At present it is here for the use of others while that process is ongoing. 
