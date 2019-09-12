#!/usr/bin/env python

"""
extract-control.py --input <infile>.gff > <outfile>.fasta

This script takes an input of a gff file (usually generated by Proka) with a
##FASTA section containing the sequence data, and outputs a FASTA file of all
the control regions found. The script is capable of finding control regions in
reverse complemented and rotated sequences, unless the feature is split across
the start/end of the sequence.

The general process is as follows:
    + Annotations for the Ile-tRNA and 12S rRNA are located.
    + If necessary, adjustments are made for reverse complemented and
        roatated sequences.
    + the appropriate section of the FASTA section is extracted and printed to
        standard output.

 !! This script has only been tested with mosquito mitogenomes, but should !!
 !! be capable of extracting control regions from any mitogenome bounded   !!
 !! by Ile-tRNA and 12S rRNA                                               !!
"""

#   Extend the system to clip any control region, given command line arguments for the annotation boundaries

# Disable circular distance?

# command line arguments, GFF parser
import argparse
from BCBio import GFF
# Type hinting
from typing import Dict, Tuple, Iterable

def get_params():
    """Returns the command line arguments."""

    parser = argparse.ArgumentParser(description="""
        This script takes an input of a gff file (usually generated by Proka) with a
        ##FASTA section containing the sequence data, and outputs a FASTA file of all
        the control regions found. The script is capable of finding control regions in
        reverse complemented and rotated sequences, unless the feature is split across
        the start/end of the sequence. Output is done through standard output.
        """.strip())
    parser.add_argument("--input", help="The file to be worked on. GFF format with a ##FASTA section.")
    parser.add_argument("--output", help="The destination of the output file.")
    parser.add_argument("--force", action="store_true", help="Overwrites the output if it already exists.")

    return parser.parse_args()

def check_flip(start_annot, end_annot) -> bool:
    """
    checks if the features have been flipped from the expected orientation (5'-B---A-3' instead of 5'-A---B-3').
    returns T/F
    """
    return (end_annot is not None
            and start_annot is not None
            and end_annot.location.start < start_annot.location.start)

def check_rev_comp(end_annot) -> bool:
    """
    Checks if the end annnotation has been flipped.
    Returns T/F
    """
    # operating under the assumption that the 12S rRNA is on the negative strand
    # so if it isn't, we assume that we're looking at the reverse complement of what we're normally looking at
    return (end_annot is not None
            and end_annot.location.strand == 1)

def parse_seq(rec, start_annot, end_annot) -> str:
    """
    Clips the relevant sequence data from sequences.

    There are four possible states for the sequence:
    - Normal orientation
    - reverse complement
    - annotations in reverse order (sequence has been rotated)
    - reverse complement rotated

    In a normal orientation, it gets sequence data from the 1 position, to the
    start of the "near" annotation, then skips all nucleotides until the end of the
    "end" annotation. These sequences are then concatenated together.

    In the reverse complement, it still grabs data from the edges, but from 1 to the
    start of the "far" annotation, and after the end of the "near" annotation. These
    annotation are only "near" and "far" in reference to the normal orientation. In
    the reverse complement, they get reversed.

    In reverse order, it gets from the end of the "far" to the start of the "near"
    sequence. Again, the annotations are only far/near in reference to the normal
    orientation of the sequence.

    In the reverse complement reverse order, it gets from the end of the "start" to
    the start of the "end." In this case, it's double reversed so the annotation names
    actually bodge their way to being correct.

    Returns a sequence object.
    """

    # return an empty sequence if one of the bounds is missing.
    if end_annot is None or start_annot is None:
        return ""

    rev_comp = False
    inner = False

    # check if it's reversed and/or flipped
    if check_rev_comp(end_annot):
        rev_comp = True
        inner = True
        # we have to set inner flag to be "flipped" because of the orientation of the sequence.

        # let me essplane:

        # Normal: 5'~tRNA------->12S>~~3'          ║ Reverse complement: 3'~~<12S<-------tRNA~5'
        # the tRNA is before the rRNA, so we take  ║ The rRNA is before the tRNA, but we still have
        # from the start, until the tRNA, then     ║ to take the outer. Because the rotation checker
        # from the rRNA to the end.                ║ only checks for the position of the rRNA vs tRNA,
        #                                          ║ it'll see the 12S before and decide that we need
        #                                          ║ to take the centre, which is false. So we pre-set
        #                                          ║ it to true to trick it into being correct by
        #                                          ║ toggling itself back to getting the edges.
        #══════════════════════════════════════════╬══════════════════════════════════════════════════
        # Rotated: 5'--->12S>~~~tRNA----3'         ║ Rotated RevComp: 3'----tRNA~~~<12S<---5'
        # The 12S rRNA is before the tRNA, so we   ║ We want the centre, but since the 12S is after
        # take the inner sequence, not the edges.  ║ the tRNA, the rotation checker will ignore it as
        # The rotation checker will see the 12S    ║ correct, and take the edges unless we set it to
        # before the tRNA and switch to taking the ║ true (which we did).
        # middle instead of the edges.             ║

        # Is this a kludgy hack born out of tech debt? Probably. Still works, tho
    if check_flip(start_annot, end_annot):
        inner = not inner # take from between instead of the tails
        # it has to be a toggle instead of an assignment because of the weirdnesss with rotated reverse complements.

    # if we're looking at the reverse complement, we need to flip the boundaries.
    # defaults
    # The "near" bound is closer to the 5' end
    near_bound = 0
    # It's the "far" bound because it's further from the 5' end.
    far_bound = len(rec)

    # Set the variables to their correct value, assuming the boundaries were found, otherwise
    # they just stay as the defaults.
    if rev_comp:
        if end_annot is not None: near_bound = end_annot.location.start
        if start_annot is not None: far_bound = start_annot.location.end
    else:
        if start_annot is not None: near_bound = start_annot.location.start
        if end_annot is not None: far_bound = end_annot.location.end

    ## DEFAULT BEHAVIOUR:
    # take from the 12S rRNA to the end, then from the beginning to the
    if not inner:
        out_seq = rec.seq[far_bound:]+rec.seq[:near_bound]

    ## MODIFIED BEHAVIOUR:
    #  Take BETWEEN the rRNA and the tRNA.
    # This is "backwards" because the expected far bound is nearer when this happens.
    else:
        out_seq = rec.seq[far_bound:near_bound]

    # if it's supposed to get reversed, reverse it.
    if rev_comp:
        out_seq = out_seq.reverse_complement()

    return str(out_seq)


def circular_distance(a: int, b: int, C: int) -> int:
    """
    Finds the shortest distance between two points along the perimeter of a circle.

    Args:
        a: a point on a circle's circumference.
        b: another point on the cicrle.
        C: the total circumference of the circle.
    Returns:
        a the shortest good.

    >>> circular_distance(2,5,10)
    3
    >>> circular_distance(12,3,15)
    6
    >>> # It even works with numbers >C or <0
    ... circular_distance(-20, 37, 10)
    3
    """
    arc = abs(a - b) % C # the distance between these in one direction -- not necessarily the shortest distance
    return min(C - arc, arc) # the arc and the complement of the arc, one of which will be shorter than the other.

def check_anchor(anchor_loc: int, seq_len: int, loc1: int, loc2: int) -> bool:
    """
    Given a location for the anchor, the total length of the sequence, and two annotations,
    Checks if the second annotation would be closer to the anchor point.
    Returns T/F

    >>> check_anchor(10,20,4,15)
    True
    >>> check_anchor(5,20,2,9)
    False
    """
    # here's how this one goes:
    # We get the circular distance of two points from the same anchor
    # Then we see if the first one is larger than the second
    return (circular_distance(loc1, anchor_loc, seq_len)
            > circular_distance(loc2, anchor_loc, seq_len))

def find_bound(seq_len: int, all_annots, anchor: int):
    """
    Takes a sequence length, a list of annotations, a tuple containing all the bounds we're looking for, and the anchor location
    for that feature.
    Finds the first annotation that matches the highest priority boundary.
    Returns an annotation or None.
    """
    best_fit = None
    for annot in all_annots:
        # if there is a best fit and it is closer than the new one, skip this annotation.
        # if there isn't a best fit, or the new annotation is closer to the anchor, reassign best_fit
        if (best_fit is None
                or not check_anchor(anchor, seq_len, int(annot.location.start), int(best_fit.location.start))):
            best_fit = annot
        if anchor is None:
            return best_fit
    # if there's an anchor, we have to get to the very end to find the best match. Otherwise, it breaks early, as above.
    return best_fit

def find_anchor(bound_start_anchor, start_anchor_annots):
    """
    Finds the first available anchor gene and returns its position
    returns the start position of the anchor or -1 if there's no anchor
    """
    for s_anch in bound_start_anchor:
        for annot in start_anchor_annots:
            if annot.qualifiers["product"][0] == s_anch[0]:
                return annot.location.start if s_anch[1] == "start" else annot.location.end # Bonus: this breaks both loops
    return None # otherwise, return the default (which is -1)

#  Dict[str, SeqFeature] ->  -> list[SeqFeature]
def get_features(product_features, products: Iterable[str]):
    return [product_features[p] for p in products if p in product_features]

def main():
    """Main CLI entry point for extract-control.py"""
    args = get_params()

    # make a new file, either forcing overwrite of the old file or not, depending on the setting.
    with open(args.input, "r") as seq_file, open(args.output, "w" if args.force else "x") as out_file:

        # highest to lowest priority, e.x. 12S is preferred to 12S (partial)
        bound_start = ("mtRNA-Ile(gat)", "mtRNA-Ile(aat)")
        bound_end = ("12S ribosomal RNA", "12S ribosomal RNA (partial)")
        # anchors are basically only for tRNAs as boundaries because tRNAs are easy to copy and move,
        # and are sometimes just totally mis-annotated when done automatically.
        # So we use a protein coding gene, preferably large so it's harder to translocate, but generally just the closest.
        # Then we find the closest tRNA to that gene and use that.
        # IF NO ANCHOR IS GIVEN OR FOUND, it just uses the first occurence of the highest priority bound marker
        bound_start_anchor = (("NADH-ubiquinone oxidoreductase chain 2","start"),("12S ribosomal RNA","end"), ("12S ribosomal RNA (partial)","end"))

        # through every sequence record in the gff
        for rec in GFF.parse(seq_file):
            # Make a dictionary called prod_features containing SeqFeature objects where the key is the product name
            # and the value is a list of features that contain the same product tags.
            prod_features: Dict[str, SeqFeature] = {f.qualifiers['product'][0]: f for f in rec.features if 'product' in f.qualifiers}

            start_anchor_annots = get_features(prod_features, bound_start_anchor)
            start_annots = get_features(prod_features, bound_start)
            end_annots = get_features(prod_features, bound_end)

            s_anchor_loc = find_anchor(bound_start_anchor, start_anchor_annots)

            # through every possible start/end bound
            start_annot = find_bound(len(rec), start_annots, s_anchor_loc)
            end_annot = end_annots[0] if len(end_annots) > 0 else None

            # get the sequence
            out_seq = parse_seq(rec, start_annot, end_annot)

            # print it out to the file
            out_file.write(f">{rec.id}_cont_reg\n{out_seq}\n")

if __name__ == '__main__':
    main()
