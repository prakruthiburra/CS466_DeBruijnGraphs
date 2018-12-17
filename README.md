# CS466_DeBruijnGraphs
Course project to analyse DeBruijn graphs for Enterobacteria phage lambda genome

PLEASE NOTE:
Parameters are currently hardcoded in the script. The code has not been modularised.

cleaning.py
This script takes a genome file with each line of the form 


1201 caaaaaacta ccgtgaaaag tcggtggatg tggcgggtta tgatgaactt gctgcttttg


It generates a file with the genome stored as a continuous string

gen_reads.py
This script generates reads given a read length, error rates and coverages.

deBruijn_size.py
This script constructs deBruijn graphs for different prefixes of the phage genome

deBruijn_genReads.py
This script generates reads from a phage genome with specified coverage, error rates and read lengths.
