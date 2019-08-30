#!/usr/bin/env python

"""
bind-gff-to-fasta.py --gff [in].gff --fasta [in2].fasta > output.gff

Takes a .gff file and a .fasta file and puts them together. Doesn't check if
the .gff already has a ##FASTA section.
"""

import argparse

def get_params():
    """Gets the command line arguments"""
    parser = argparse.ArgumentParser(description="""
        Takes a .gff file and a .fasta file and puts them together. Doesn't check if
        the .gff already has a ##FASTA section. Output is done through standard
        output.
        """.strip())
    parser.add_argument("--gff", help="The gff file to be worked on.")
    parser.add_argument("--fasta", help="The FASTA sequence file to be worked on. Make sure it matches the annotation file.")

    return parser.parse_args()

def main():
    """
    this just prints all the contents of the .gff, then prints the ##FASTA tag, then all the fasta sequences.
    it should probably go directly into a file instead of relying on standardoutput, but whatever, it works.
    """
    args = get_params()

    with open(args.gff) as gff_file:
        for i in gff_file:
            print(i, end='')

    print("##FASTA")

    with open(args.fasta) as fasta_file:
        for i in fasta_file:
            print(i, end='')

if __name__ == '__main__':
    main()
