# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 09:41:50 2022

@author: Victoria Gayle Buswell

The code below was used to perform the abba baba anaysis on a dataframe of allele frequencies created from a vcf file. 
Original file paths have been removed and lists emptied ready for use. The last three loops are the ones that perform all test. You need to comment and uncomment to use them. 
"""
import pandas as pd
import numpy as np
import scipy
from scipy import stats

#sample lists, these are lists of the names of the samples that represent each

high_car = []#insert full list of low introgression A. m. carnica (0.99 membership to correct cluster)
high_car_a = []#insert half of the low introgression A. m. carnica (0.99 membership to correct cluster)
high_car_b = []#insert half of the low introgression A. m. carnica (0.99 membership to correct cluster)

south_west = [] #insert list of all the samples you want to test. For example all south west samples

high_amm =[]#insert full list of low introgression A. m. mellifera (0.99 membership to correct cluster)
high_amm_a = []#insert half of the low introgression A. m. mellifera (0.99 membership to correct cluster)
high_amm_b = []#insert half of the low introgression A. m. mellifera (0.99 membership to correct cluster)

high_lig = []#insert full list of low introgression A. m. ligustica (0.99 membership to correct cluster)
high_lig_a =  []#insert half of the low introgression A. m. ligustica (0.99 membership to correct cluster)
high_lig_b = [ ]#insert half of the low introgression A. m. ligustica (0.99 membership to correct cluster)



# choose your Dstat p1 p2 and p3 populations from the lists above.
# You need an additional populations to be used in the f stat for p3a and p3b.
#this is why you split the list at the top, so for example, high_car_a and high_car_b can act as p3a and p3b:



def d_stat(p1, p2, p3):
    '''Returns D statistic based on three samples.
    Postive D indicates introgression between p2 and p3.
    Negative D introgresison between p1 and p3.
    Does not indicate the proportion of introgression'''
    ABBA = (1- p1)*p2*p3
    BABA = p1*(1 - p2)*p3
    num = sum(ABBA) - sum(BABA)
    denom = sum(ABBA) + sum (BABA)
    d = num / denom
    return d

def f_stat(P1, P2, P3a, P3b):
    abba_num = (1-P1)*P2*P3a
    baba_num = P1*(1-P2)*P3a
    abba_denom = (1-P1)*P3b*P3a
    baba_denom = P1*(1-P3b)*P3a
    numer = (sum(abba_num)) - (sum(baba_num))
    denom = (sum(abba_denom)) - (sum(baba_denom))
    f = numer / denom
    return f

def block_jackknife_per_snp(frequ_table, p1, p2, p3, out):
    '''Returns d statistic block jackknife procedure.
    Blocks are based on number of SNPs per block.'''
    with open(out, 'w') as d_out:
        row_start = 0
        row_end = 525
        while row_end > 0:
            df = pd.read_csv(frequ_table, delim_whitespace=True)
            row_total = df.shape[0]
            temp_df = df.drop(df.index[row_start:row_end])
            d_answer = d_stat(temp_df[p1],temp_df[p2],temp_df[p3])
            d_out.write('SNPs '
                        +'\t'
                        +str(len(temp_df))
                        +'\t'
                        +str(d_answer)
                        +'\n')
            row_start += 525
            row_end += 525
            if row_end > row_total:
                temp_df = df.drop(df.index[row_start:row_total])
                d_answer = d_stat(temp_df[p1],temp_df[p2],temp_df[p3])
                d_out.write('SNPs'
                            +'\t'
                            +str(len(temp_df))
                            +'\t'
                            +str(d_answer)
                            +'\n')
                break


def significant_d(allele_freq ,d_stat_results, d_significant, p1,p2,p3):
    '''This function will take the output from block jackknife functions.
    Then calculate standard and z scores, and p-values for the d statistics'''
    with open(d_significant, 'w') as sig_out:
        df_d = pd.read_csv(d_stat_results, sep='\t', header =None)
        print(df_d)
        df_afs = pd.read_csv(allele_freq, delim_whitespace=True)
        print(df_afs)
        d = d_stat(df_afs[p1],df_afs[p2],df_afs[p3])
        print(d)
        d_sd = df_d[2].std()
        variance = np.var(df_d[2])
        print(d_sd)
        number_of_rows = df_d.shape[0]
        var_blocks = variance*number_of_rows
        norm_dis_err = np.sqrt(var_blocks)
        print(number_of_rows)
        print('Approx normally distributed SE  = '
              +str(norm_dis_err))
        #root = np.sqrt(number_of_rows)
        #d_err = d_sd/root
        #print(d_err)
        d_z = d / norm_dis_err
        print('Z - score = '
              + str(d_z))
        #D_p <- 2*pnorm(-abs(D_Z)) <--- this R function in python? to get p-value of D
        #pnorm is cumulative density function (cdf.)
        #abs in r computes the absolute value.
        #scipy.stats.norm.cdf() is the equivelent of pnorm in python
        p = 2*scipy.stats.norm.cdf(-abs(d_z))
        #these both come out the same, both ways of calculating a p-value.
        p = scipy.stats.norm.sf(abs(d_z))*2
        print('P value is: '
              +str(p))
        sig_out.write('P1:P2:P3'+'\t'+'Overall D'+'\t'+'Standard d'+'\t'+'Normally distributed standard error'+'\t'+'Z score'+'\t'+'P value'+'\n'
                      +P1+':'+P2+':'+P3+'\t'+str(d)+'\t'+str(d_sd)+'\t'+str(norm_dis_err)+'\t'+str(d_z)+'\t'+str(p))


def admixture(allele_frequencies_in, p1, p2, p3a, p3b, results_out):
    df_afs = pd.read_csv(allele_frequencies_in, delim_whitespace=True)
    print(df_afs)
    f_result = f_stat(df_afs[p1], df_afs[p2], df_afs[p3a], df_afs[p3b])
    print(f_result)
    with open(results_out, 'w') as out:
        out.write(str(p1)+
                  ';'+
                  str(p2)+
                  ';'+
                  str(p3a)+
                  '\t'+
                  str(p3b)+
                  '\t'+
                  str(f_result))
        
#below are three sets of loops to perform block jack knifing and significance test as well as f statistics
#comment and uncomment as needed:


# This loop and functions perform significance testing for trio MEL;SW;CAR. WARNING, results in lots of files, place them somewhere sensisble!:

# for p_one in high_amm:
#     for p_two in south_west:
#         for p_three in high_car:
#             P1 = p_one
#             P2 = p_two
#             P3 = p_three
#             file_suffix = f'trios_p1_{P1}_p2_{P2}_p3_{P3}_20Blocks_525SNPs_amm_sw_car'
#             AF = '\\path\\to\\your\\colony_allelefreq_dataframe'
#             OUT_SNP = '\\pah\\you\\want\\the\\blockjackknife\\out\\files_%s'  % file_suffix
#             IN_SNP ='\\path\\to\\your\\blockjackknife\\results\\files_%s'  % file_suffix
#             OUT_SIG = '\\path\\to\\where\\you\\want\\the\\d_significant_stats_%s.txt'  % file_suffix
#             block_jackknife_per_snp(AF, P1, P2, P3, OUT_SNP) #this is the blockjackknife function as above
#             significant_d(AF,IN_SNP, OUT_SIG, P1, P2, P3)  #this is performed on the blockjack knife out files to estimate significance for each trio
            
# This loop and functions perform significance testing for trio MEL;SW;LIG, WARNING, results in lots of files, place them somewhere sensisble!:

# for p_one in high_amm:
#     for p_two in south_west:
#         for p_three in high_lig:
#             P1 = p_one
#             P2 = p_two
#             P3 = p_three
#             file_suffix = f'trios_p1_{P1}_p2_{P2}_p3_{P3}_20Blocks_525SNPs_amm_sw_lig'
#             AF = ''\\path\\to\\your\\colony_allelefreq_dataframe.txt'
#             OUT_SNP = '\\pah\\you\\want\\the\\blockjackknife\\out\\'  % file_suffix
#             IN_SNP ='\\path\\to\\your\\blockjackknife\\results\\files_%s.txt'  % file_suffix
#             OUT_SIG = '\\path\\to\\where\\you\\want\\the\\d_significant_stats_%s.txt'  % file_suffix
#             block_jackknife_per_snp(AF, P1, P2, P3, OUT_SNP)
#             significant_d(AF,IN_SNP, OUT_SIG, P1, P2, P3)  
  
#    This loop results in f, proportion of admixture testing: 

# for p_one in high_amm:
#     for p_two in south_west:
#         for p_three in high_lig_a:
#             for p_three_b in high_lig_b:
#                 file_suffix = f'trios_p1_{p_one}_p2_{p_two}_p3_{p_three}_f_amm_sw_lig_lig'
#                 AF = '\\path\\to\\your\\colony_allelefreq_dataframe.txt'
#                 out = '\\path\\to\\where\\you\\want\\the\\f_result_%s' % file_suffix
#                 admixture(AF, p_one, p_two, p_three, p_three_b, out)
                

            
            
            
            
            
            
            
            
            