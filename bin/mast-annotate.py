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


def get_params():
    """Gets command line arguments. Returns them."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--mast", help="The motif.xml file to be worked on.")

    return parser.parse_args()

class Seq: #sequence
    """
    Stores some information about the sequence: name and length.
    """

    def __init__(self, n, l):
        self.name = n
        self.length = l


class Mot: #motif
    """
    Stores motif data. name, altname if it has one, and the length of the motif
    """
    def __init__(self, n, a, l):
        self.name = n
        self.altn = a
        self.length = l

def main():
    """He does it"""
    args = get_params()

    seqs = {} # sequences
    mot = {}  # motifs
    ann = {}  # annotations

    # MAST uses a "reverse complement" tag, so no = + strand, yes = - strand
    strand = {
        "n": "+",
        "y": "-"
    }

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
            seqs[len(seqs)] = Seq(seq_tag["name"], int(seq_tag["length"]))

            # then, go through every hit under each sequence to find the actual motif locations.
            for hit_tag in seq_tag.find_all("hit"):
                cur_motif = mot[int(hit_tag["idx"])]

                # [sequence 1 name]      MEME Suite                    nucleotide_motif
                # [position]            [position + motif-length]      .
                # [strand]              .                              ID=[ID];Name=[motifName]

                # from this, we construct the annotations.
                ann[len(ann)] = (seq_tag["name"]+                                   # sequence name
                                 "\tMEME Suite"+                                                 # source
                                 "\tnucleotide_motif\t"+                                         # type
                                 str(int(hit_tag["pos"])+1)+"\t"+                                # start
                                 str(int(hit_tag["pos"])+cur_motif.length)+                      # end
                                 "\t."+                                                          # score
                                 "\t"+strand[hit_tag["rc"]]+                                     # strand
                                 "\t."+                                                          # frame
                                 "\tNote=p-value:"+hit_tag["pvalue"]+";Name="+cur_motif.altn+cur_motif.name)  # note & other

    # This could be slightly shorter if
    print("##gff-version 3")
    for i in seqs:
        print("##sequence-region "+seqs[i].name+" 1 "+str(seqs[i].length))
    for i in ann:
        print(ann[i])

# import-safety
if __name__ == '__main__':
    main()
