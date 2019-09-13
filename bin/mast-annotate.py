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

    parser.add_argument("--input", help="The motif.xml file to be worked on.")
    parser.add_argument("--output", help="The output destination (usually a .gff)")
    parser.add_argument("--force", action="store_true", help="Overwrites the output if it already exists.")

    return parser.parse_args()

def get_xml_data(xml_path):
    """
    Gets the XML data from the file and closes it nicely, returning the output.
    """
    with open(xml_path) as mast_file:
        return BeautifulSoup(mast_file, 'lxml')

def get_seq_info(xml_path: str):
    """

    """
    seqs = [] # sequences
    mot = []  # motifs

    # crack the input file open and look at it 👀
    xml = get_xml_data(xml_path)

    # grab attributes from all the motif tags
    for mot_tag in xml.find_all("motif"):
        # make a new object in the table
        mot.append(mot_tag)

    # grab all the sequence tags and get their names and lengths.
    for seq_tag in xml.find_all("sequence"):
        seq = SeqRecord("A"*(int(seq_tag["length"])),seq_tag["name"])
        seqs.append(seq)

        # then, go through every hit under each sequence to find the actual motif locations.
        for hit_tag in seq_tag.find_all("hit"):
            cur_motif = mot[int(hit_tag["idx"])]

            # from this, we construct the annotations.
            qualifiers = {
                "source": "MEME Suite",
                "Note": "p-value:"+hit_tag["pvalue"],
                # If the alt name is empty, just the name. Otherwise, altn+" "+name
                "Name": " ".join(filter(None, [cur_motif.attrs.get("alt", ""), cur_motif["id"]]))
            }

            # build the sequence feature object
            seq.features.append(
                SeqFeature(
                    FeatureLocation(
                        int(hit_tag["pos"])-1, # MAST indexes at 1, biopython indexes at 0
                        int(hit_tag["pos"])+int(cur_motif["length"])-1
                    ),
                    type="nucleotide_motif",
                    strand=-1 if hit_tag["rc"] == "y" else 1,
                    qualifiers=qualifiers
                )
            )

    return seqs

def main():
    """Main CLI entry point for mast-annotate.py"""
    args = get_params()

    # open the outfile and start writing the sequence list to the file
    with open(args.output, "w" if args.force else "x") as out_file:
        GFF.write(get_seq_info(args.input), out_file)

if __name__ == '__main__':
    main()
