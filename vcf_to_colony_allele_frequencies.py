# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 09:04:26 2021

@author: Victoria Gayle Buswell
"""
import sys
import io
import os
import pandas as pd
import numpy as np

def get_AD(column):
    all_sites = []
    for x in column:
        z= x.split(':')
        y = z[1].split(',')
        #print(z[2])
        try:
            ref = (float(y[0]) / float(z[2]))
        except ValueError:
            ref = 0
        except ZeroDivisionError:
            ref = 0
        try:
            alt = (float(y[1]) / float(z[2]))
        except ValueError:
            alt = 0
        except ZeroDivisionError:
            alt = 0
        all_sites.append(ref)
        all_sites.append(alt)
    return(all_sites)     
        

def read_vcf(path,af_out,SNP_num):
    out_file = open(af_out, 'w')
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    vcf = pd.read_csv(
        io.StringIO(''.join(lines)), delim_whitespace=True
    ).rename(columns={'#CHROM': 'CHROM'})
    print(vcf)
    SNP_ref = SNP_num - 1
    new_df = vcf.filter(['CHROM','POS','REF','ALT'])
    df_2 = pd.concat([new_df]*2, ignore_index = True)
    df_2['Ref/Alt']= df_2.iloc[:SNP_ref, 4] = 'ref' # row up to the number of SNPs in final file and column 4 
    df_2.iloc[SNP_num:,4] = 'alt' # row from SNPs +1 then column 4 is alt
    for column in vcf.iloc[:,9:]: #starting at the first sample column (position9) in the vcf and go to the last column 
        allele_frequencies= get_AD(vcf[column]) # get allele frequencies from the AD field
        af =pd.Series(allele_frequencies) # make the the af a series
        df_2[column]= af.values # add it the the data frame
    df_2.to_csv(out_file, sep='\t', index =True)# print the new allele frequency data frame out. 
    
    

#read_vcf('E:\\path\\to\\your\\VCF_file.vcf','E:\\path\\to\\your\\allele_freq_out_file.txt', 15000)

file_in=str(sys.argv[1])
file_out=str(sys.argv[2])
snp_num=int(sys.argv[3])
read_vcf(file_in, file_out, snp_num)
