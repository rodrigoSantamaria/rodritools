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

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re
import os


fastaPath=None
outPath=None
wigPath=None
threshold=None
window=150
step=50
above=True
#%%
def getValue(arg, val):
    if(val[0]=='-'):
        print("ERROR: a value expected for argument ", arg, " (use -h for help).")
        sys.exit()
    return val
    
def printHelp():
    print('\ngcWindows version 0.1')
    print('\nDetermines CG content from sequences [and correlates with occupancy]')
    print('\nUsage: ')
    print('\tpython gcWindows.py -fasta file [-wig file] [-out file] [options]')
    print('\t-fasta file:\tFasta file with sequences to extract fragments.')
    print('\t-wig file:\tWig file with nucleosome occupancy.')
    print('\t-out file:\tName of the bed file were fragments will be saved (bed format)')
    print('\t\t(Default name is inputFastaFile.tsv)')
    print('\nAdditional options:')
    print('\t-t threshold:\tIf set, only windows with GC percentage above wich fragments are extracted (default 75).')
    print('\t-w window:\tFragment size (default 150).')
    print('\t-s step:\tStep among fragments (default 50).')
    print('\t-l:\t Extracts fragments *below* the GC threshold (default above).')
    
    print('')
     
    print('\t-h:\t\tPrint this help.')

    print('\nRodrigo Santamaría, et al. rodri@usal.es, Dpto, de Informática y Automática, Universidad de Salamanca.')
    print('\t\tDistributed under GPL3.0 License')
    print('')
#%%
if __name__ == "__main__":
#    print("ARGV", sys.argv)
    
    if len(sys.argv)==1:
        print('\t\tgcWindows version 0.1')
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
        if sys.argv[i]=="-wig": 
            wigPath=getValue(sys.argv[i], sys.argv[i+1])
        if sys.argv[i]=="-out":
            outPath=getValue(sys.argv[i], sys.argv[i+1])
        if sys.argv[i]=="-t":
             threshold=int(getValue(sys.argv[i], sys.argv[i+1]))
        if sys.argv[i]=="-w":
             window=int(getValue(sys.argv[i], sys.argv[i+1]))
        if sys.argv[i]=="-s":
             step=int(getValue(sys.argv[i], sys.argv[i+1]))
        if sys.argv[i]=="-l":
             above=False
           
    if outPath==None:
        outPath=fastaPath+".tsv"

    #fastaPath="/home/rodri/data/annotations/sp/Spombe.fasta"
    #wigPath="/home/rodri/data/nucleosomas/Dhrp3_MNase/Dhrp3_Rep_1_depth_wl_trimmed_PE.wig"
    print("READING INPUT...")
    seqs=bioio.readFasta(fastaPath)
    if(wigPath!=None):
        wig=bioio.readWig(wigPath)
    
    #%% Operation
    print("EXTRACTING FRAGMENTS...")
    print("[Fragment size="+str(window)+", threshold>="+str(threshold)+"%]")
    
    fragments={}
    gcContent=[]
    wigValues=[]
    for k in seqs.keys():
        print(k)
        s=seqs[k]
        for i in range(0,len(s)-window,step):
            content=100*(bioseq.gcContent(s[i:(i+window)])/float(window))
            #print(int(content))
            if(threshold==None or (above and content>threshold) or ((not above) and content<=threshold)):
                fragments[k+":"+str(i)+"-"+str(i+window)]={"chromosome":k, "start":i,"end":i+window,"score":int(content)}
                if(wigPath!=None):
                    value=np.mean(wig[k][i:i+window])
                    wigValues.append(value)
                    fragments[k+":"+str(i)+"-"+str(i+window)]["name"]=value
                    #print(fragments[k+":"+str(i)+"-"+str(i+window)])
                gcContent.append(content)
                
    if(len(wigValues)==0 and len(gcContent)!=len(wigValues)):
        print("The size of wig and fasta files is different")
        
    #%% Saving
    print(str(len(fragments))+" fragments extracted.")
    bioio.writeBedFromDict(outPath, fragments)
    
if(False):    
    #%% Graphs
    #bin occupancy by gcContent
    bins=[]
    for i in range(99):
        bins.append([])
    for i in range(len(gcContent)):
            bins[int(gcContent[i])].append(wigValues[i])
    mbins=np.zeros(100)
    sbins=np.zeros(100)
    for i in range(99):
        mbins[i]=(np.mean(bins[i]))
        sbins[i]=(np.std(bins[i]))
    #%%
    zeroVals=[]
    wigValues2=[]
    gcContent2=[]
    for i in range(len(wigValues)):
        if(wigValues[i]==0):
            zeroVals.append(i)
        else:
            gcContent2.append(gcContent[i])
            wigValues2.append(wigValues[i])
    #%%
    #b,m=np.polyfit(gcContent2, np.log(wigValues2), 1)
    z=np.polynomial.polynomial.polyfit(gcContent2, np.log(wigValues2),1)
    p=np.poly1d(z)
    #p=np.poly1d()
    plt.plot(gcContent2, np.log(wigValues2), 'bo', markersize=1)   
    plt.plot(gcContent2, p(gcContent2))
    plt.xlabel("GC Content")
    plt.ylabel("Average occupancy (log)")
    #%%
    plt.scatter(gcContent, np.log(wigValues), marker='.')
    #%%
    sns.regplot(gcContent, np.log(wigValues), markersize=1)
    #%%
    plt.plot(range(100), mbins, label="GC content")
    plt.fill_between(range(100), mbins-sbins, mbins+sbins)
    plt.xlabel("GC Content")
    plt.ylabel("Average occupancy")
    