# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 10:57:47 2021

@author: Victoria Buswell
"""
import sys

def imiss_filter(filein, missing, fileout):
    filein = open(filein, "r")
    fileout = open(fileout, 'w')
    row = filein.readlines()[1:]
    print(row)
    for line in row:
        item = line.split()
        print(item[4])
        if float(item[4]) > missing:
            fileout.write(f'{item[0]}\n')
    
    
filein=str(sys.argv[1])
missing= float(sys.argv[2])
fileout=str(sys.argv[3])
imiss_filter(filein, missing, fileout)
