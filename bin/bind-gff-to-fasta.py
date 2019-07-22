#!/usr/bin/env python3 
# 
# .gff specialization
# Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'
#    seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
#    source - name of the program that generated this feature, or the data source (database or project name)
#    feature - feature type name, e.g. Gene, Variation, Similarity
#    start - Start position of the feature, with sequence numbering starting at 1.
#    end - End position of the feature, with sequence numbering starting at 1.
#    score - A floating point value.
#    strand - defined as + (forward) or - (reverse).
#    frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
#    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
# 
# [sequence 1 name]	MEME Suite	nucleotide_motif	[position]	[position + motif-length]	.	[strand]	.	ID=[ID];Name=[motifName]
# 
# 


## Is it bad?
## Yes

## Does it do what I want it to?
## Also yes


import argparse
import re

def getParams():
	parser = argparse.ArgumentParser()
	parser.add_argument("--gff", help = "The file to be worked on")
	parser.add_argument("--fasta", help = "The file to be worked on")

	return parser.parse_args()

def main():

	## Tag-specific search functions
	# These get passed as arguments sometimes because if it allows this kind of weirdness, I'm down.
	# putting them in here is easier than passing around all the variables as arguments

	# if it's a <sequence/> line, parse it into a seq object
	
	args = getParams()

	gff_file = open(args.gff)

	for i in gff_file:
		line = i[:-1] # clip the newline off
		print(line)

	gff_file.close()

	print("##FASTA")
	
	fasta_file = open(args.fasta)

	for i in fasta_file:
		line = i[:-1]
		print(line)

	fasta_file.close()

# import-safety
if __name__ == '__main__':
	main()