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
filteredPath=None
delete=False
rev=False
t=95
#%%
def getValue(arg, val):
    if(val[0]=='-'):
        print("ERROR: a value expected for argument ", arg, " (use -h for help).")
        sys.exit()
    return val
    
def printHelp():
    print('\nsimilarSequences version 0.1')
    print('\nRetrieves sequences with high homology')
    print('\nUsage: ')
    print('\tpython similarSequences.py -fasta file [OPTIONS]')
    print('\t-fasta file:\tFasta file with sequences to search for (returns a file.tsv  with the sequence pairs above the homology threshold.')
    print('\nOptions:')
    print('\t-out file:\tFasta file with sequences to search for.')
    print('\t-t XX:\tIdentity threshold (0-100). Only sequence pairs above this threshold will be reported (default 95, no decimals allowed).')
    print('\t-d:\tDeletes the first occurrence of similar sequences above t, saving the filtered sequences on a file with name "input file-filter-t.fasta"')
    print('\t-r:\tDetermines similarity also with reverse complement sequences"')
    
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
        if sys.argv[i]=="-d": 
            delete=True
        if sys.argv[i]=="-r": 
            rev=True
        if sys.argv[i]=="-t":
             t=int(getValue(sys.argv[i], sys.argv[i+1]))
             if(t>100 or t<0):
                 print("Threshold must be a number between 0 and 100")
                 exit()
                
    if delete:
        filterPath=fastaPath.replace(".fasta", "-filter"+str(t)+".fasta")
    
    if outPath==None:
        outPath=fastaPath.replace(".fasta","")+".tsv"

    print("Reading sequences...")
    #fastaPath="/home/rodri/data/nucleosomas/Nuc_SP.fasta"
    #fastaPath="/home/rodri/workspace/python/ibfg/rodritools/test.fasta"
    
    seqs=bioio.readFasta(fastaPath)
    
    
    #%% Operation
    print("Comparing sequences (it may take several minutes)....")
    print(rev, t)
    aln=bioseq.identity(seqs,t, revComp=rev)
    #%% Save results
    print("Saving results to "+outPath)
    f=open(outPath,"w")
    f.write("seq1\tseq2\tidentity\n")
    for k in aln.keys():    
        #print(k)
        kk=k.split(",")
        f.write(kk[0]+"\t"+kk[1]+"\t"+str(round(100*aln[k]))+"\n")
        #print(kk[0]+"\t"+kk[1]+"\t"+str(round(100*aln[k]))+"\n")
    f.close()  
    #%%
    seqs2=seqs.copy()
    if delete:
        for k in aln.keys():
            k0=k.split(",")[0]
            if k0 in seqs2.keys():
                del seqs2[k0]
        bioio.writeFasta(seqs2, filterPath)
    