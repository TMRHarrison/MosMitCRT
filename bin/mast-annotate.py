#!/usr/bin/env python
#
# .gff specification
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
# [sequence 1 name] MEME Suite  nucleotide_motif    [position]  [position + motif-length]   .   [strand]    .   ID=[ID];Name=[motifName]
#

"""
mast-annotate.py --mast [mast].xml > [output].gff

Takes a MAST output .xml file and makes a gff annotation file out of it.
"""


import argparse
from bs4 import BeautifulSoup

from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def get_params():
    """Gets command line arguments. Returns them."""
    parser = argparse.ArgumentParser(description="""
        Takes a MAST output .xml file and makes a gff annotation file out of it.
        """.strip())

    parser.add_argument("--mast", help="The motif.xml file to be worked on.")
    parser.add_argument("--out", help="The output destination (usually a .gff)")

    return parser.parse_args()

class Mot: #motif
    """
    Stores motif data. name, altname if it has one, and the length of the motif
    """
    def __init__(self, n, a, l):
        self.name = n
        self.altn = a
        self.length = l

def main():
    """Main CLI entry point for mast-annotate.py"""
    args = get_params()

    seqs = [] # sequences
    mot = {}  # motifs

    # MAST uses a "reverse complement" tag, so no = + strand, yes = - strand
    strand = {
        "n": -1,
        "y": 1
    }

    # crack the input file open and look at it 👀
    with open(args.mast) as mast_file:
        xml = BeautifulSoup(mast_file, 'lxml')

        # grab attributes from all the motif tags
        for mot_tag in xml.find_all("motif"):
            alt_name = ""
            #if there is an alt name, add it
            if mot_tag.has_attr("alt"):
                alt_name = mot_tag["alt"]+" "

            # make a new object in the table
            mot[len(mot)] = Mot(mot_tag["id"], alt_name, int(mot_tag["length"]))

        # grab all the sequence tags and get their names and lengths.
        for seq_tag in xml.find_all("sequence"):
            seqs.append(SeqRecord("A"*(int(seq_tag["length"])),seq_tag["name"]))

            # then, go through every hit under each sequence to find the actual motif locations.
            for hit_tag in seq_tag.find_all("hit"):
                cur_motif = mot[int(hit_tag["idx"])]

                # from this, we construct the annotations.
                qualifiers = {
                    "source": "MEME Suite",
                    "Note": "p-value:"+hit_tag["pvalue"],
                    "Name": cur_motif.altn+cur_motif.name
                }

                # build the sequence feature object
                seqs[len(seqs)-1].features.append(
                    SeqFeature(
                        FeatureLocation(
                            int(hit_tag["pos"])-1, # MAST indexes at 1, biopython indexes at 0
                            int(hit_tag["pos"])+cur_motif.length-1
                        ),
                        type="nucleotide_motif",
                        strand=strand[hit_tag["rc"]],
                        qualifiers=qualifiers
                    )
                )

    # open the outfile and start writing the sequence list to the file
    with open(args.out, "w") as out_file:
        GFF.write(seqs, out_file)

if __name__ == '__main__':
    main()
