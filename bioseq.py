# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 09:58:56 2017

@author: rodri
"""

#%%---------------------- INVERSA COMPLEMENTARIA ----------------------
#Retorna la cadena (inversa) complementaria a seq
def complement(seq):
    comp=""
    complementos={"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}#tabla de conversión
    for s in seq.upper()[::-1]:                 #recorrer en sentido inverso
        comp=comp+complementos[s]
    return comp

#%%--------------------- INVERSA ------------------------------
def inverted(seq):
    inv=""
    for s in seq.upper()[::-1]:
        inv=inv+s
    return inv

#------------------------------- ALIGNMENTS -------------------------------
#%% Get the consensus sequence for a dictionary of sequences
#The consensus letters are shown in lowercase if their frequency is in the range [lower,upper)
#and in upper case if are above the upper parameter. If below the 'lower' parameter noconsensus character is given
def consensus(seqs, lower=.5, upper=.99, noconsensus="-", weights={'A':1,'C':1,'G':1,'T':1,}):
    match=[]
    lseq=len(list(seqs.values())[0])
    for j in range(lseq):
        freq=dict(A=0,G=0,T=0,C=0,N=0)
        for s in seqs.values():
            freq[s[j]]+=1*weights[s[j]]
       # print(freq)
       # add=max(freq.iteritems(), key=operator.itemgetter(1)) #not working in python 3
        add=max(freq, key=(lambda k:freq[k]))
        maxf=max(freq.values())/len(seqs.values())
        if(maxf>upper):
            match.append(add.upper())
        elif(maxf>lower):
            match.append(add.lower())
        else:
            match.append(noconsensus)
    return ''.join(match) #return reduce(lambda a,b: '{0}{1}'.format(a,b), match)  #HIC SUNT DRACONES: más info sobre las funciones lambda: http://www.secnetix.de/olli/Python/lambda_functions.hawk

#%% 4) Compute frequencies 
# seqs - an array of nucleotide sequences of the same length (uppper case)
# returns a dictionary with nucleotides as keys and arrays of frequencies for 
#  the key, with a length equal to the length of each sequence in seqs
def frequency(seqs, k=1):
    import numpy as np
    
    sl=len(seqs[0])
    fr={}
    for k in ["A","C","T","G", "N"]:
        fr[k]=np.zeros(sl)
    for i in range(sl): #for each nucleotide
        for s in seqs:      #for each sequence
            fr[s[i]][i]+=1
    for k in fr.keys():
        for i in range(len(fr[k])):
            fr[k][i]=(float)(fr[k][i])/len(seqs)
    return fr

#motifs must be a dictionary with sequences as values
#Returns the number of nucleotides that differ from the consensus
def score(motifs):
    s=0
    con=consensus(motifs)
    for m in motifs.values():
        for j in range(len(m)):
            if(con[j].upper()!=m[j].upper()):
                s+=1
    return s
    
#score(alnSeqs)
#cs=consensus(remseqs)  
#%%      
"""Given an dict of aligned sequences (aln), returns the pairwise identities
of the sequences, without counting gaps (-) or not-known nucleotides (N) on
any of the sequences. 
"""
def identity0(aln):
    ret={}
    for i in range(len(aln)):
        k1=list(aln.keys())[i]
        s1=list(aln.values())[i]
        for j in range(i+1, len(aln)):
            k2=list(aln.keys())[j]
            s2=list(aln.values())[j]
            s=0
            slen=0
            for p in range(len(s2)):
                if( (s1[p].upper()!="N" and s1[p]!="-") and (s2[p].upper()!="N" and s2[p]!="-") ):
                    slen=slen+1
                    if(s1[p].upper()==s2[p].upper()):
                        s+=1
            #print(k1+", "+k2+"\t"+str(round(100*s/slen,2)))
            if(slen>0):
                ret[k1+","+k2]=s/slen
            else:
                ret[k1+","+k2]=0
    return ret

#.09s  por seq (comparando contra 75K seqs, normal e inversa)
def identity(aln, t=95, revComp=True):
    ret={}
    import time
    keys=list(aln.keys())
    values=[x.upper() for x in list(aln.values())]
    maxDif=len(values[0])*(1-t/100.0)
    alnc={}
    for k in aln.keys():
        alnc[k]=complement(aln[k])
    t0=time.time()
    for i in range(len(keys)):
        #print(i)
        #t00=time.time()#10e-5
        k1=keys[i]
        s1=values[i]
        for j in range(i+1, len(aln)):
            k2=keys[j]
            s2=values[j]
            d=0
            for p in range(len(s2)):
                if(s1[p]!=s2[p]):
                    d+=1
                    if(d>maxDif):
                        break
            if(d<maxDif):
                ret[k1+","+k2]=(len(s2)-d)/len(s2)
            if(revComp and d>=maxDif):
                d=0
                s2r=alnc[k2]
                for p in range(len(s2r)):
                    if(s1[p]!=s2r[p]):
                        d+=1
                        if(d>maxDif):
                            break
                if(d<maxDif):
                    ret[k1+","+k2]=(len(s2r)-d)/len(s2r)
        #print(time.time()-t00)
    print("time in ", time.time()-t0)
    return ret
#tal=identity(seqs, revComp=True)
#alnSeqs={'s1':"AAATTT", 's2':"TTTAAA", 's3':"TTTTTT", 's4':"AAAAAA"}
#tal=identity(alnSeqs, revComp=True,t=80)
    


#%%------------------------ CODON MANIPULATION -------------------------
#Returns a dictionary where keys are the 3-nucleotide sequence and value the corresponding aminoacid
def codonTable():    
    bases = ['T', 'C', 'A', 'G']
    codons = [a+b+c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))
    return codon_table

#Returns a dictionary where keys are the aminoacid and values the list of all 3-nuc corresponding sequences
def inverseCodonTable():
    import numpy as np
    ct=codonTable()
    ict={}
    for aa in np.unique(list(ct.values())):
        ict[aa]=list(filter(lambda k:ct[k]==aa, list(ct.keys())))
    return ict

#Takes a nucleotide sequence and alters it so it contains the same aminoacid
#sequence but mutates randomly across the possible codons
def randomCodons(seq):
    import random
    ct=codonTable()
    ict=inverseCodonTable()
    rseq=""
    for i in range(0,len(seq),3):
        cod=seq[i:i+3]
        rcod=random.sample(ict[ct[cod]],1)[0]
        rseq+=rcod
    return rseq

#%% ---------------------- COMPARISON ----------------------------
#Returns the homology between two same-lenght nucleotide sequences s1 and s2 as 
#a dictionary:
#   aa - amino acid homology (considering a 0-frame start)
#   codon - codon homology (considering a 0-frame start)
#   nt - nucleotide homology
def coincidence(s1,s2,aspects=["codon", "nt", "aa"]):
    co={}
    
    if("codon" in aspects):
        co["codon"]=sum([1 if s1[i:i+3]==s2[i:i+3] else 0 for i in range(0,len(s1),3)])/(len(s1)/3.)
    if("nt" in aspects):
        ls1=len(s1)
        co["nt"]=sum([1.0 if s1[i]==s2[i] else 0.0 for i in range(ls1)])/ls1
    if("aa" in aspects):
        ct=codonTable()
        ict=inverseCodonTable()
        co["aa"]=sum([1 if (s1[i:i+3] in ict[ct[s2[i:i+3]]])==True else 0 for i in range(0,len(s1),3)])/(len(s1)/3.)
    return co
    
def coincidenceNT(s1,s2):
    ls1=len(s1)
    return sum([1.0 if s1[i]==s2[i] else 0.0 for i in range(ls1)])/ls1
    
#------------------- IUPAC AMBIGUOUS CODES -------------------------
#%% Given an IUPAC sequence, it unfolds it on all its possible nucleotide combinations
def unfold(seq):
   from Bio import Seq
   import itertools
   d = Seq.IUPAC.IUPACData.ambiguous_dna_values
   r = []
   for i in itertools.product(*[d[j] for j in seq]):
      r.append("".join(i))
   return r
#%%
"""
Given an expression with degenerate nucleotides, it provides a regular expression
for searches
seq - sequence with nucleotides (including degenerate ones (B-Y)
returns a regular expression where each nucleotide has been substituted by
        all the nucleotides it comprises (e.g. AWR --> A[ATW][AGR])
"""
def nuc2exp(seq):
    #mapper according to http://opisthokonta.net/?p=549
    nmap={"B":"[CGTBSKY]",
         "D":"[AGTDRWK]",
         "H":"[ACTHMYW]",
         "K":"[GTK]",
         "M":"[ACM]",
         "N":"[ACGTBDHKMNRSVWY]",
         "R":"[AGR]",
         "S":"[CGS]",
         "V":"[ACGVMSR]",
         "W":"[ATW]",
         "Y":"[CTY]",
         "A":"A",
         "C":"C",
         "G":"G",
         "T":"T",
         }
    exp=""
    for s in seq:
        exp+=nmap[s]
    return exp



#%% ------------------------- MUTATIONS ---------------------
# Dada una secuencia (word), un número de mutaciones puntuales (num_mismatches)
# y una cadena con los posibles componentes (letters), devuelve una variable 
# generator con todas las posibles mutaciones. 
# OJO: el resultado hay queconvertirlo a list para su porterior uso
#Código extraido de: 
#http://stackoverflow.com/questions/11679855/introducing-mutations-in-a-dna-string-in-python
def mutations(word, num_mismatches, letters="ACGT"):
    import itertools
    for locs in itertools.combinations(range(len(word)), num_mismatches):
        this_word = [[char] for char in word]
        for loc in locs:
            orig_char = word[loc]
            this_word[loc] = [l for l in letters if l != orig_char]
        for poss in itertools.product(*this_word):
            yield ''.join(poss)
            
            
# ------------------------- MUTACIONES (II) ---------------------
# Como la función anterior, pero ahora retorna una lista con las mutaciones 
#posibles con num_mismathes SNPs *o menos*
def mutationsEqualOrLess(word, num_mismatches, letters="ACGT"):
   matches=set()
   for dd in range(num_mismatches,-1,-1): 
       matches.update(list(mutations(word, dd, letters)))
   return matches

# Pruebas, cuántas mutaciones de 2/3 SNPs hay en un 9-mer?
# kmer='CACAGTAGGCG'
# mu=list(mutations(kmer, 2))
# print len(mu)
# mu=list(mutations(kmer, 3))
# print len(mu)
# mu=mutationsEqualOrLess(kmer, 3)
# len(mu)

# ------------------------- MUTACIONES (III) --------------------- S2E5
# Retorna como lista todos los posibles {k}-mers en la str {seq} y 
#todas sus mutaciones en hasta 2 SNPs
def allMutations(seq, k=9, d=2):
    kmers=set()
    for i in range(len(seq)-k+1):
        kmer=seq[i:i+k]
        kmers.update(mutationsEqualOrLess(kmer, d))
    return kmers


#%% Mutates a nucleotide sequence seq in k positions.
# Only A, T,C, G allowed as characters in the sequence.
# If repeatPos is True, some of these positions may be mutated several times (default False). 
#If repeatNuc is True, a nucleotide may be mutated into itself (default False)
# w sets the weights or probability (as number adding up to 100) of choosing each
# nucleotide.
def mutate(seq, k=1, repeatPos=False, repeatNuc=False, w={"A":25, "C":25, "G":25, "T":25}):
    import random 
    pos=list(range(len(seq)))
    nuc=["A","T","C","G"]
    ww=[]
    for n in nuc:
        ww.append(w[n])
    rseq=list(seq.upper())
    for i in range(k):
        #select a random position
        rpos=random.choice(pos)
        if not repeatPos:
            pos.remove(rpos)
        #select a random (maybe weighted) nucleotide
        ant_nuc=rseq[rpos]
        rseq[rpos]=random.choices(nuc, weights=ww, k=1)[0]
        if not repeatNuc:
            while(ant_nuc==rseq[rpos]):
                rseq[rpos]=random.choices(nuc, weights=ww, k=1)[0]
                
            
    return "".join(rseq)        
#tal=mutate("TTT", k=3, repeatNuc=True)   
#print(tal)     
    
#%% Sampling based on respecting {k}-mer frequencies
def sample(seq,k=1, replace=False):
    import numpy as np
    options=[]
    for i in range(len(seq)-k+1):
        options.append(seq[i:i+k])
    return "".join(np.random.choice(options, int(len(seq)/k),replace=replace))
#mut=sample("AATTCCGG", 2)
#%%
def atContent(seq):
    return seq.count("A")+seq.count("T")

def gcContent(seq):
    return seq.count("G")+seq.count("C")

#%%
"""
Given a {seq} (str or array of str), returns the number of occurrences of each
dinucleotide included in the list.
If raw=False, returns its total sum instead
"""
def dnContent(seq, dn=[], raw=True):
    cont={}
    for k in dn:
        cont[k]=0
    if(type(seq)==str):
        seq=[seq]
    for s in seq:
        for k in cont:
            cont[k]+=s.count(k)
    if(not raw):
        return sum([x for x in cont.values()])
    else:
        return cont
#print(dnContent("ATGTGTGGCGCCGCACGCTA", ["GC","CG","CC","GG"]))
    
#%%
"""
Given a {seq} (str or array of str), returns the number of occurrences of each
{k}-mer. If {freq} is True (default), it returns instead the % of each k-mer
"""
def content(seq, k=1, freq=True):
    comb=list(allMutations("A"*k, k=k, d=k))
    cont={}
    for k in comb:
        cont[k]=0
    if(type(seq)==str):
        seq=[seq]
    for s in seq:
        for k in cont:
            cont[k]+=s.count(k)
    if not freq:
        return cont
    else:
        fr={}
        totalCont=sum(cont.values())
        for k in cont:
            fr[k]=round(100*cont[k]/totalCont,2)
        return fr
#%%
def nuc2aa(seq):
    ct=codonTable()
    aa=""
    for i in range(0,len(seq),3):
        aa=aa+ct[seq[i:i+3]]
    return(aa)
#%%
def shuffle(seq):
    import random
    rseq=[*seq]
    random.shuffle(rseq)
    return("".join(rseq))


""" Returns a dicts with {num} sequences of length {size} taken from random positions on
{genome}, which is a dict as the one returned by bioio.readFasta.
Sequences that contain Ns are discarded.
"""
def randomSequencesByLocation(genome, num=10000, size=50):
    import random
    supergenome=""
    for k in genome:
        supergenome=supergenome+genome[k]
    gseqs={}
    for i in range(num):
        gseqs[i]="NN"
        while("N" in gseqs[i]):
            ri=random.randint(0, len(supergenome))
            gseqs[i]=supergenome[ri:ri+size]
    return gseqs

""" Returns a dict with {num} sequences of length {size} formed by random nucleotides, 
chosen with even probabilities or from the probabilities in a given
{genome}, which is a dict as the one returned by bioio.readFasta.
"""
def randomSequencesByContent(genome=None, num=10000, size=50, k=1):
    import bioseq
    import random
    
    weights={}
    if(genome!=None):
        supergenome=""
        for track in genome:
            supergenome=supergenome+genome[track]
        weights=bioseq.content(supergenome, k)
    else:
        comb=list(allMutations("A"*k, k=k, d=k))
        w=1./len(comb)
        for key in comb:
            weights[key]=w
        #weights={"A":0.25,"C":.25,"G":.25,"T":.25}
    rseqs={}
    for i in range(num):
        rseqs[i]="".join(random.choices(list(weights.keys()), weights=list(weights.values()), k=int(size/k)))
    return rseqs
# tal=randomSequencesByContent(sp_genome, num=1, size=10000,k=2)
# supergenome=""
# for track in sp_genome:
#     supergenome=supergenome+sp_genome[track]
# con=bioseq.content(supergenome,2)



#%% Methods:
"""
{nucseq} must be a dictionary with sequences as values.
returns a dictionary with the periodicity per dinuclotide of each {dinucs} pair
combination
if {genome}!=None, a normalization is done based on random sequences
"""
def periodicity(nucseq, dinucs=["AA", "TT"], maxdist=51, genome=None):
    import time
    import numpy as np
    dro={}
    cont=0
    keys=[]
    normalize=False
    if(genome!=None):
        normalize=True
        rseq=randomSequencesByLocation(genome=genome, num=len(nucseq), size=150)

    for dn1 in dinucs:
        for dn2 in dinucs:
            #print(dn1,"-",dn2)
            print(cont)
            #test
            if complement("".join([dn1,dn2])) not in keys:
                keys.append("".join([dn1,dn2]))
                print(keys[-1])
                ro=np.zeros(maxdist)
                for seq in nucseq.values():
                    #seq=seq[75:]
                    t0=time.process_time()
                    for i in range(len(seq)-len(dn1)+1):
                        if(seq[i:i+len(dn1)]==dn1):
                            for j in range(i+1,min(i+maxdist,len(seq)-len(dn1)+1)):
                                if(seq[j:j+len(dn2)]==dn2):
                                    ro[j-i]+=1
                #random control
                if(normalize):
                    rro=np.zeros(maxdist)
                    for seq in rseq.values():
                        t0=time.process_time()
                        for i in range(len(seq)-len(dn1)+1):
                            if(seq[i:i+len(dn1)]==dn1):
                                for j in range(i+1,min(i+maxdist,len(seq)-len(dn1)+1)):
                                    if(seq[j:j+len(dn2)]==dn2):
                                        rro[j-i]+=1
                    #normalization
                    dro["-".join([dn1,dn2])]=[ro[i]/max(1,rro[i]) for i in range(len(ro))]
                else:
                    dro["-".join([dn1,dn2])]=ro
                cont+=1
    return dro
#%% Determine minima and maxima in a dict of numeric arrays
def localMinMax(ro, prominence=None, d=None):
    import scipy.signal
    peaks={}
    for k in ro:
        peaks[k]={}
        peaks[k]["position"]=scipy.signal.find_peaks(ro[k], distance=d, prominence=prominence)[0]
        peaks[k]["prominence"]=scipy.signal.peak_prominences(ro[k], peaks[k]["position"], wlen=10)
        peaks[k]["width"]=scipy.signal.peak_widths(ro[k], peaks[k]["position"], prominence_data=peaks[k]["prominence"], wlen=20)
    return peaks     


#%% Visualization
def plotPeriodicity(dro, titleadd="", name="", lim=None, ylines=None):
    import matplotlib.pyplot as plt
    for k in dro.keys():
        nro=dro[k]
        plt.plot(range(2,len(nro)), nro[2:], "-")
        if(titleadd!=""):
            plt.title(titleadd+"-"+k)
        else:
            plt.title(k) 
        if(lim!=None and len(lim)==2):
            plt.ylim(lim)
        if(ylines!=None):
            for yl in ylines[k]["position"]:
                if(yl>2):
                    plt.plot([yl,yl], [min(nro[2:]), max(nro[2:])], "--", color="grey",linewidth=1)
                    plt.text(yl+0.2,min(nro[2:]), yl)
        plt.xlabel("Distance")
        plt.ylabel("Occurrences")
        if(name!=""):
            plt.savefig(name+"-"+k+".png")
        else:
            plt.show()
        plt.clf()
        