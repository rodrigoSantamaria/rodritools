#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Codon synonym shuffling
Given a coding sequence, select synonymous codons and shuffle them across positions
That is, if sequence is 

TCT TCC CCG TCG
Ser Ser Pro Ser

Get a random permutation of synonym codons such as:

TCC TCG CCG TCT
Ser Ser Pro Ser
    
@author: rodri
"""

#Methods for I/O
#import os
#os.chdir("/home/rodri/workspace/python/ibfg/biotools")
import bioio
import bioseq

import random
#%%-------------------------------------------------------------------------
#Shuffling method (requires a dict of seqs as returned by readFasta)
"""
if strict=True, it would not replace each position by the same codon at 
that position (notice that it does not prevent to replace it by another
identical codon at another position)
"""
def shuffleSynonymCodons(seqs, strict=True):
    #codon tables
    ict=bioseq.inverseCodonTable()
    ct=bioseq.codonTable()
    
    #initialize dict for storing syn codons and positions
    scod={}
    for k in ict.keys():
        scod[k]=[]
        
    #1st phase: collect synonym codons and positions
    for k in seqs.keys():
        seq=seqs[k]
        for i in range(0, len(seq), 3):
            cod=seq[i:i+3]
            scod[ct[cod]].append((i,cod))
            
    
    #print(scod)
    #2nd phase: shuffle synonym codons
    nseqs={}
    for k in seqs.keys():
        #print(k)
        seq=seqs[k]
        nseq=""
        for i in range(0, len(seq), 3):
            cod=seq[i:i+3]
            aa=ct[cod]
            ri=(i,cod)
            
            print(i,"-------------------------------------------")
            print("opciones: ", scod[aa], "\noriginal", ri)
            if(len(scod[aa])>1):
                while(ri[0]==i):
                    ri=random.sample(scod[aa], k=1)[0] 
                scod[aa].remove(ri)
            else:
                ri=scod[aa][0]
                scod[aa]=[]
            print("elección: ",ri)
            ncod=ri[1]
            nseq=nseq+ncod
        nseqs[k]=nseq
        
    return nseqs
            
#%Usage: enter to readFasta the path to your ORF file
s=bioio.readFasta("/home/rodri/data/ura4/ura4orf.fasta")
sseq=shuffleSynonymCodons(s)

#%%Write to another fasta:
#bioio.writeFasta(sseq, "/home/rodri/data/ura4/ura4orf-shuffledSynCod.fasta")
s1=s["ura4_pombe"]
s2=sseq["ura4_pombe"]
#% Checks: AT content, nucleotide content, amino acid sequence
print(bioseq.atContent(s["ura4_pombe"]))
print(bioseq.atContent(sseq["ura4_pombe"]))


print(bioseq.nucContent(s["ura4_pombe"]))
print(bioseq.nucContent(sseq["ura4_pombe"]))

print(bioseq.nuc2aa(s["ura4_pombe"]))
print(bioseq.nuc2aa(sseq["ura4_pombe"]))
aa1=bioseq.nuc2aa(s["ura4_pombe"])
aa2=bioseq.nuc2aa(sseq["ura4_pombe"])
print(aa1==aa2)
