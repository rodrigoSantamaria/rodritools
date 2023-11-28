#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 12:46:35 2022
BWT searches based on bwsearchMethods.py

@author: rodri
"""


#%%
import sys

import bioio
import bioseq
import bwsearchMethods

import numpy as np
import re
import os
import pickle


fastaPath=None
outPath=None
indexPath=None
saveIndex=False
index=None
seq=None
rep=None
repPath=None
d=0
#%%
def getValue(arg, val):
    if(val[0]=='-'):
        print("ERROR: a value expected for argument ", arg, " (use -h for help).")
        sys.exit()
    return val
    
def printHelp():
    print('\nbwsearch version 0.1')
    print('\nProvides fast inexact sequence searches')
    print('\nUsage: ')
    print('\tpython bwsearch.py -fasta file -seq sequence [-out file] [options]')
    print('\t-fasta file:\tFasta file with sequences to search for.')
    print('\t-seq sequence:\tSequence to search for.')
    print('\t-rep sequence:\tSequence to replace for occurrences (by default no replacement done).')
    print('\t-out file:\tName of the bed file were occurrences will be saved. The score column contains the sequence (bed format)')
    print('\t\t\tIf -rep is especified, a fasta file with the same name of the bed file is also generated')
    print('\t\t(if not provided, default name is fasta-seq-d.bed)')
    print('\nAdditional options:')
    print('\t-d mutations:\tNumber of mutations allowed over the sequence (default 0, max allowed 2).')
    print('\t-i:\tIf provided, the generated index is saved with name fastaFile.bwt')
    print('\t-ii inputIndex:\tIf provided, the index is not generated but read from this file.')
    #print('\t-s step:\tStep among fragments (default 50).')
    #print('\t-l:\t Extracts fragments *below* the GC threshold (default above).')
    
    print('')
    
    print('\t-h:\t\tPrint this help.')

    print('\nRodrigo Santamaría, et al. rodri@usal.es, Dpto, de Informática y Automática, Universidad de Salamanca.')
    print('\t\tDistributed under GPL3.0 License')
    print('')
#%%
if __name__ == "__main__":
#    print("ARGV", sys.argv)
    
    if len(sys.argv)==1:
        print('\t\tbwsearch version 0.1')
        print('\t\tFor usage info, please try: python bwsearch.py -h')
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
            repPath=outPath.replace(".*$","")+"-replaced.fasta"
        if sys.argv[i]=="-seq":
            seq=getValue(sys.argv[i], sys.argv[i+1])
        if sys.argv[i]=="-rep":
            rep=getValue(sys.argv[i], sys.argv[i+1])
        if sys.argv[i]=="-d":
             d=int(getValue(sys.argv[i], sys.argv[i+1]))
        if sys.argv[i]=="-i":
             saveIndex=True
        if sys.argv[i]=="-ii":
             indexPath=getValue(sys.argv[i], sys.argv[i+1])
             f=open(indexPath,"rb")
             index=pickle.load(f)
            
    if(saveIndex):
        indexPath=fastaPath.replace(".fasta", ".bwt")
  
    if outPath==None:
        outPath=fastaPath.replace(".fasta","")+"-"+seq+"-"+str(d)+".bed"
        repPath=fastaPath.replace(".fasta","")+"-"+seq+"-"+str(d)+"-replaced.fasta"

    seqs=bioio.readFasta(fastaPath)
    
    if(d>0 and len(seq)/d<5):
        print("The ratio between the length of the searched pattern (-seq) and the number of allowed mutations (-d) must be kept above 4 to make fast searches")
    
    #%% Operation
    result=bwsearchMethods.searchLocal(seqs, seq, d, indexPath=indexPath, index=index)
    result=result["points"]
   
    #%% Save results
    print("Saving results to "+outPath)
    f=open(outPath,"w")
    for k in result.keys():
        for x in result[k]:
            f.write(str(k)+"\t"+str(x)+"\t"+str(x+len(seq))+"\t.\t"+seqs[k][x:x+len(seq)]+"\t.\n")
    f.close()  
    #%% If replacement
    seqs2=seqs.copy()
    if rep != None:
        for k in result.keys():
            temp=list(seqs2[k])
            for x in result[k]:
                temp[x:(x+len(seq))]=rep
            seqs2[k]="".join(temp)
    bioio.writeFasta(seqs2, repPath)