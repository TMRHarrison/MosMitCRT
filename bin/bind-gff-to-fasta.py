#!/usr/bin/env python

## Is it bad?
## Yes

## Does it do what I want it to?
## Also yes

import argparse

def getParams():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", help = "The gff file to be worked on.")
    parser.add_argument("--fasta", help = "The FASTA sequence file to be worked on. Make sure it matches the annotation file.")

    return parser.parse_args()

# this just prints all the contents of the .gff, then prints the ##FASTA tag, then all the fasta sequences.
# it should probably go directly into a file instead of relying on standardoutput, but whatever, it works.
def main():

    args = getParams()

    with open(args.gff) as gff_file:
        for i in gff_file:
            line = i[:-1] # clip the newline off
            print(line)

    print("##FASTA")

    with open(args.fasta) as fasta_file:
        for i in fasta_file:
            line = i[:-1]
            print(line)

# import-safety
if __name__ == '__main__':
    main()