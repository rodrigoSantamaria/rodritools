#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 10:05:35 2022

@author: rodri
"""

#%%
import sys

import bioio
import bioseq

import numpy as np
import re
import os


fastaPath=None
outPath=None
decimals=2
#%%
def getValue(arg, val):
    if(val[0]=='-'):
        print("ERROR: a value expected for argument ", arg, " (use -h for help).")
        sys.exit()
    return val
    
def printHelp():
    print('\ngcContent version 0.1')
    print('\nUsage: ')
    print('\tpython gcContent.py -fasta file [-out file]')
    print('\t-fasta file:\tFasta file with sequences to compute their GC content.')
    print('\t-out folder:\tName of the file were GC content will be saved.')
    print('\t\tEach line contains the GC content (as a value in [0-1]) of the corresponding sequence.')
    print('\t\t(Default name will be inputFastaFile-content.tsv.')
    print('\t-d decimals:\tPrecision for the contents (default 2).')
    print('')
     
    print('\t-h:\t\tPrint this help.')

    print('\nRodrigo Santamaría, et al. rodri@usal.es, Dpto, de Informática y Automática, Universidad de Salamanca.')
    print('\t\tDistributed under GPL3.0 License')
    print('')
#%%
if __name__ == "__main__":
#    print("ARGV", sys.argv)
    
    if len(sys.argv)==1:
        print('\t\tgcContent version 0.1')
        print('\t\tFor usage info, please try: python gcContent.py -h')
        sys.exit()

    #Checking for required parameters
    if ('-h' in sys.argv)==False:
        if ('-fasta' in sys.argv)==False:
            print("ERROR: an input file must be provided (use -h for help).")
            sys.exit()
            
        
    #input retrieval
    for i in range(1, len(sys.argv)):
        
        if sys.argv[i]=='-h' or len(sys.argv)==1:
            printHelp()
            sys.exit()

        #input params
        if sys.argv[i]=="-fasta": 
            fastaPath=getValue(sys.argv[i], sys.argv[i+1])
        if sys.argv[i]=="-out":
            outPath=getValue(sys.argv[i], sys.argv[i+1])
            
    if outPath==None:
        outPath=fastaPath.replace(".fasta","")+"-content.tsv"

    #fastaPath="/home/rodri/data/ndr/ndrSP/version1/NDRs_genes_Spombe.fasta"
    seqs=bioio.readFasta(fastaPath)
    #outFile="/home/rodri/data/ndr/ndrSP/version1/NDRs_genes_Spombe.fasta"
    
    #%% Operation
    contents={}
    for k in seqs.keys():
        s=seqs[k]
        contents[k]=round(bioseq.gcContent(s)/len(s),decimals)
    #%%
    f=open(outPath,"w")
    f.write("seq\tGC_content\n")
    for k in contents.keys():
        f.write(k+"\t"+str(contents[k])+"\n")
    f.close()
    