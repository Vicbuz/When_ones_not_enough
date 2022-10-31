# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 12:44:09 2021

@author: Victoria Gayle Buswell
"""


import csv
import pandas as pd
import numpy as np
def iPLEX_to_ped(in_file, new_ped):
    ped = open(new_ped,"w")
    with open(in_file) as file:
        iPLEX = pd.read_csv(in_file)
        iPLEX = iPLEX.fillna('0 0') #replace the empty genotypes with the number 0 for ped
        # take all monomorphic sites and add the other base, as all markers nust be biallelic
        iPLEX_G = iPLEX.replace('G', 'G G').replace('GA','G A').replace('GC','G C').replace('GT','G T')
        iPLEX_A = iPLEX_G.replace('A','A A').replace('AG', 'A G').replace('AT','A T').replace('AC','A C')
        iPLEX_T = iPLEX_A.replace('T', 'T T').replace('TA','T A').replace('TG','T G').replace('TC','T C')
        iPLEX_C= iPLEX_T.replace('C','C C').replace('CA','C A').replace('CG','C G').replace('CT','C T')
        iPLEX_final = iPLEX_C
        print(iPLEX_final) # iPLEX FINAL has all the biallelc changes

        pivot = iPLEX_final.pivot(index="Sample Id", columns="Assay Id")# create a wide not long dataframe
        print(pivot)
        pivot2 = pivot.drop(['Description'], axis=1)
        listofzeros = [0] * 62 #this repesents the 62 samples in the CSV but would be any number!
        #62 samples in csv
        listofsamples= #insert a list of sample names for example listofsamples = [1,6,11,16,21,26,31,36,41,42,43,216,221,226,231,236]
        pivot2.insert(loc=0, column='paternal', value=listofzeros)
        pivot2.insert(loc=0, column='maternal', value=listofzeros)
        pivot2.insert(loc=0, column='sex', value=listofzeros)
        pivot2.insert(loc=0, column='not_used', value=listofzeros)
        pivot2.insert(loc=0, column='sample', value=listofsamples)
        
        pivot2.to_csv(new_ped, index=True, header=True)

        
    ped.close()

iPLEX_to_ped('/path/to/your/csv.csv', '/out/path/to/your/new/ped.txt')
