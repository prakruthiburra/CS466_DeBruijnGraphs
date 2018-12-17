# CS466_DeBruijnGraphs
Course project to analyse DeBruijn graphs for Enterobacteria phage lambda genome

PLEASE NOTE:
Parameters are currently hardcoded in the script. The code has not been modularised.
Code to generate de Bruijn graphs has been taken from http://nbviewer.jupyter.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_deBruijn.ipynb 
Genome source - https://www.ncbi.nlm.nih.gov/nuccore/215104 


1. cleaning.py
This script takes a genome file with each line of the form 


1201 caaaaaacta ccgtgaaaag tcggtggatg tggcgggtta tgatgaactt gctgcttttg


It generates a file with the genome stored as a continuous string

2. gen_reads.py
This script generates reads given a read length, error rates and coverages.

3. deBruijn_size.py
This script constructs deBruijn graphs for different prefixes of the phage genome

4. deBruijn_genReads.py
This script generates reads from a phage genome with specified coverage, error rates and read lengths.
