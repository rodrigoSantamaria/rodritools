#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 08:56:36 2018

This program converts files from baMM format to meme format
BAMM
.ihbcp -> format with 1 'paragraph' by nucleotide in the motif.
First line on the paragraph is single nucleotide frequencies (4 numbers, 
supposed in order ACGT)
Second line seems th have dinucleotide freqs (16 numbers, order unknown)
Third line has 64 numbers so it should be trinucleotide freqs.

I'd ignore the 2nd and 3rd line as meme only works with single nt.

There's 1 file per motif here
MEME
.meme has a header which we can easily replicate
then it has a 'background letter freqs' paragraph that i'd mimick with .25 for 
each nt

Then there comes the MOTIFs paragraphs with
1) motif name 
2) leter prob matrix info (i'd mimick)
3) motif matrix - correspondent to the bamm first lines

@author: rodri
"""
import sys
#%%
def bamm2meme(bammFolder="", memeFile=""):
    import re
    if(memeFile==""):
        memeFile=bammFolder+".meme"
    
    from os import listdir
    from os.path import isfile, join
    onlyfiles = [f for f in listdir(bammFolder) if isfile(join(bammFolder, f))]
    ihbcps=sorted([bammFolder+"/"+f for f in onlyfiles if f.endswith(".ihbcp")])
    print(ihbcps)
    f=open(memeFile, "w")
    f.write("MEME version 4\n\n")
    f.write("ALPHABET=ACGT\n\n")
    f.write("Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
    
    
    for file in ihbcps:
        motname=re.sub("\\..*$", "", re.sub("^.*/", "", file))
        f.write("MOTIF "+motname+"\n")
        f.write("unknown, retrieved from BaMM file "+file+"\n")
        f.write("letter-probability matrix:")
        #f.write("letter-probability matrix: alength=4, w=")
        fm=open(file)
        motmat=[]
        for line in fm.readlines():
            if(line.count(" ")==4):
                motmat.append(line)
        #f.write(str(len(motmat))+"\n")
        f.writelines(motmat)
        f.write("\n")
    f.close()    
    return "meme sucessfuly built at"+memeFile

#%%

def printHelp():
    print('\nbamm2meme version 0.1')
#    print('For help information for each function, try:\npython remaster.py <function> -h'
    print('\nUsage: ')
#    print('\tpython remaster.py -bed|-fasta|-freqs file [-seqs file] [-ofreqs file] [-rseqs file] [-img folder] [...]'
    print('\tpython bamm2meme.py bammFolder [memeFile]')
    print('\n\t--- Arguments ---')
    print('\t-bammFolder:\tPath to a baMM folder as returned by this tool (must contain an .ihbcp file per motif)')
    print('\tmemeFile:\tPath were the converted meme format (v 4) file will be saved. By default it will be stored in the same path as de BaMM folder, with its same name and .meme extension')
   
    print('\t-h:\t\tPrint this help.')

    print('\nRodrigo Santamaría, et al. rodri@usal.es, Dpto, de Informática y Automática, Universidad de Salamanca.')
    print('\t\tDistributed under GPL3.0 License')
    print('')
#%%
if __name__ == "__main__":
#    print("ARGV", sys.argv)
    
    if len(sys.argv)<=1:
        print("ERROR: bammFolder not provided")
        printHelp()
        sys.exit()
    memeFile=""
    if(len(sys.argv)==3):
        memeFile=sys.argv[2]

    bamm2meme(sys.argv[1], memeFile)
