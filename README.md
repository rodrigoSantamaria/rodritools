# rodritools
Python methods for genomic sequence manipulations
The package consists on several python scripts with methods ready to use. Methods are commented within each script, describing arguments and output.
A minimal description of what you can find on each script can be found below

## Install
You can just download the whole package to your `destination_folder`. On python, you can just do something like:
```
import os
os.chdir(destination_folder)
import bioio
...
bioio.readFasta(...)
```

## bioio
Methods for **bio**informatics **i**nput and **o**utput
There are methods for reading fasta, wig, bed, gff, go or goa files (as well as some others). They usually require the path to the file as input and return a dict or list with the information
Analogously, there are methods for writing some of these formats.

There are also some conversion methods to convert bed to fasta or fasta to bed.

## bioseq
Methods for basic **seq**uence manipulation, such as:
* Reverse complement sequences
* Consensus sequences
* k-mer frequency occurrences
* codon table manipulations
* sequence comparison (nucleotide, codon and aminoacid homology)
* IUPAC sequence unfolding
* Generate mutated sequences
* GC, dinucleotide and k-mer content
* Sequence randomization (either as random positions from reference genomes, or synthetic sequences based on k-mer frequencies)
* Sequence periodicity

## biosearch
Methods for **search**es over sequences and other data
* Search for genes in GFFs, and obtain their sequences, given the reference genome
* Search for descriptions in GO
* Search for a given motif on sequences (this last one is ready to be used in a terminal with `python biosearch.py`)

## bwsearch
Given a fasta file, and a sequence, it returns all the positions where the sequence is in the fasta, allowing up to d mutations.
The method uses BWT to speed up the searches.
Can be used from a terminal with `python bwsearch.py` (use `h` for options)
You can also use the methods, which are available at `bsearchMethods.py`

## dnContent

## dnWindows

## frecuencias

