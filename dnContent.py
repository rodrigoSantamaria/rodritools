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
dnlist=["GC","CC","GG","CG"]
#%%
def getValue(arg, val):
    if(val[0]=='-'):
        print("ERROR: a value expected for argument ", arg, " (use -h for help).")
        sys.exit()
    return val
    
def printHelp():
    print('\ndnContent version 0.1')
    print('\nUsage: ')
    print('\tpython dnContent.py -fasta file [-out file] [options]')
    print('\t-fasta file:\tFasta file with sequences to compute their dinucleotide content.')
    print('\t-out folder:\tName of the file were content will be saved.')
    print('\t\tEach line contains the content (as a value in [0-1]) of the corresponding sequence.')
    print('\t\t(Default name will be inputFastaFile-content.tsv.')
    print('\t-dn is a comma separated list of dinucleotides to count (default GC,CC,GG,CG)')
    print('\t-d decimals:\tPrecision for the contents (default 2).')
    print('')
     
    print('\t-h:\t\tPrint this help.')

    print('\nRodrigo Santamaría, et al. rodri@usal.es, Dpto, de Informática y Automática, Universidad de Salamanca.')
    print('\t\tDistributed under GPL3.0 License')
    print('')
#%%
if __name__ == "__main__":
    
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
        if sys.argv[i]=="-d":
            decimals=int(getValue(sys.argv[i], sys.argv[i+1]))
        if sys.argv[i]=="-dn":
            dnlist=getValue(sys.argv[i], sys.argv[i+1]).split(",")
            print(dnlist)
            
    if outPath==None:
        outPath=fastaPath.replace(".fasta","")+"-content.tsv"

    seqs=bioio.readFasta(fastaPath)
    
    #%% Operation
    contents={}
    for k in seqs.keys():
        s=seqs[k]
        contents[k]={}
        counts=bioseq.dnContent(s, dnlist)
        for dn in dnlist:
            contents[k][dn]=round(counts[dn]/(len(s)-1),decimals)
        contents[k]["TOTAL"]=round(sum(list(counts.values()))/(len(s)-1),decimals)
        
    #%%
    f=open(outPath,"w")
    f.write("seq\t"+"\t".join(dnlist)+"\tsum\n")
    for k in contents.keys():
        tw="\t".join([str(x) for x in list(contents[k].values())])
        f.write(k+"\t"+tw+"\n")
    f.close()
    