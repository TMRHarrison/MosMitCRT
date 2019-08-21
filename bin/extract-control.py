#!/usr/bin/env python

# TODO:
#   Make the sequences print to files instead of being held in memory?
#   Extend the system to clip any control region, given command line arguments for the annotation boundaries

## Naming:
# variable_names
# functionNames
# Classes

# only have to do this so python 2 is compatible. When can py2 die please
import sys
if (sys.version_info[0] == 2):
    from string import maketrans

import argparse
import re

def getParams():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help = "The file to be worked on. GFF format with a ##FASTA section.")

    return parser.parse_args()

class SeqInfo:
    start_anno_5_prime = -1
    start_anno_3_prime = -1
    end_anno_5_prime = -1
    end_anno_3_prime = -1
    clipped_before = ""
    clipped_after = ""
    inner_flag = False
    rev_comp = False

# check if the features has been flipped from the expected orientation (5'-B---A-3' instead of 5'-A---B-3')
def checkFlip(cur_seq):
    ## safety for sequences when they're both found, and the END sequence is after
    if (cur_seq.start_anno_5_prime != -1
    and cur_seq.end_anno_5_prime != -1
    and cur_seq.end_anno_5_prime < cur_seq.start_anno_5_prime):
        cur_seq.inner_flag = (not cur_seq.inner_flag) # take from between instead of the tails
        # it has to be a toggle instead of an assignment because of weirdnesss with rotated reverse complements.

# biopython can also do this, but I feel like it's easier to do it this way if I just need the one thing.
# makes a translation dictionary in case of reverse compliments
nt_all  = "AGCTURYSWKMBVDH" # everything else will get ignored -> '-' and 'N' get "complemented" to themselves
nt_comp = "TCGAAYRSWMKVBHD"
# actually complements reverses the sequence
def reverseComp(seq):
    # So, because some part the python community refuses to upgrade to Python3, here we go.
    if (sys.version_info[0] == 3):
        nt_comp_dict = seq.maketrans(nt_all+nt_all.lower(), nt_comp+nt_comp.lower())
    else: #elif caveman:
        nt_comp_dict = maketrans(nt_all+nt_all.lower(), nt_comp+nt_comp.lower())

    return seq.translate(nt_comp_dict)[::-1]

def main():

    args = getParams()

    # The annotation you want the cut to START at
    bound_start_pattern = re.compile("product=mtRNA-Ile\(...\)") ## escaped regex
    # The annotation the cut ENDS at
    bound_end_pattern = re.compile("product=12S ribosomal RNA")

    # dictionary of sequences indexed by name
    seqs = {}

    # current base position being parsed
    current_base = -1

    # name of the sequence being clipped
    clipping_seq = None


    with open(args.input,'r') as seq_file:
        for line in seq_file:
            line = line[:-1] # clip the newline off
            words = line.split(None, 8) # split by whitespace

            # find the species name if you don't have one
            if words[0] == "##sequence-region":
                # Uses the sequence name as a reference so we don't have to loop
                # through the whole set every time we encounter something
                seqs[words[1]] = SeqInfo()

            # find the indices to start and stop at
            if words[0] in seqs:
                cur_seq = seqs[words[0]]
                # Sections of the gff format (to make this less opaque):
                # [0] sequence ID
                # [1] source (usually a program name)
                # [2] feature type
                # [3] start of feature
                # [4] end of feature
                # [5] score (???)
                # [6] + or - strand
                # [7] phase (reading frame) relative to the start of the annotation
                # [8] attributes: other information

                ## only use the first tRNA-Ile
                if (cur_seq.start_anno_5_prime < 0
                and bound_start_pattern.search(words[8])):
                    cur_seq.start_anno_5_prime = int(words[3]) - 1
                    cur_seq.start_anno_3_prime = int(words[4])
                    checkFlip(cur_seq)
                if bound_end_pattern.search(words[8]):
                    cur_seq.end_anno_5_prime = int(words[3]) - 1
                    cur_seq.end_anno_3_prime = int(words[4])
                    # operating under the assumption that the 12S rRNA is on the negative strand
                    # so if it isn't, we assume that we're looking at the reverse complement of what we're normally looking at
                    if words[6] == "+":
                        cur_seq.rev_comp = True
                        # we have to set inner flag to be "flipped" because of the orientation of the sequence.
                        cur_seq.inner_flag = True

                        # let me essplane:
                        #
                        # Normal: 5'~tRNA------->12S>~~3'          | Reverse complement: 3'~~<12S<-------tRNA~5'
                        # the tRNA is before the rRNA, so we take  | The rRNA is before the tRNA, but we still have
                        # from the start, until the tRNA, then     | to take the outer. Because the rotation checker
                        # from the rRNA to the end.                | only checks for the position of the rRNA vs tRNA,
                        #                                          | it'll see the 12S before and decide that we need
                        #                                          | to take the centre, which is false. So we pre-set
                        #                                          | it to true to trick it into being correct by
                        #                                          | toggling itself back to getting the edges.
                        #------------------------------------------+--------------------------------------------------
                        # Rotated: 5'--->12S>~~~tRNA----3'         | Rotated RevComp: 3'----tRNA~~~<12S<---5'
                        # The 12S rRNA is before the tRNA, so we   | We want the centre, but since the 12S is after
                        # take the inner sequence, not the edges.  | the tRNA, the rotation checker will ignore it as
                        # The rotation checker will see the 12S    | correct, and take the edges unless we set it to
                        # before the tRNA and switch to taking the | true (which we did).
                        # middle instead of the edges.             |


                        # Is this a kludgy hack born out of tech debt? Probably. Still works, tho

                    checkFlip(cur_seq)

            # when you reach a different sequence, start processing it
            if words[0][:1] == ">":
                current_base = 0
                clipping_seq = words[0][1:]

            # if the line doesn't start with >, and clipCounter >= 0
            # then you start
            elif current_base >= 0:
                # if we're looking at the reverse complement, we need to flip the boundaries.
                if seqs[clipping_seq].rev_comp:
                    near_bound = seqs[clipping_seq].end_anno_5_prime
                    far_bound = seqs[clipping_seq].start_anno_3_prime
                else:
                    near_bound = seqs[clipping_seq].start_anno_5_prime
                    far_bound = seqs[clipping_seq].end_anno_3_prime

                ## DEFAULT BEHAVIOUR:
                #  take from before the tRNA and then after the rRNA
                # clip from the start until whenever the tRNA begins
                if not seqs[clipping_seq].inner_flag:
                    seqs[clipping_seq].clipped_before += words[0][
                    :max(0,near_bound - current_base)]

                    # Clip after getting to the end of the rRNA
                    seqs[clipping_seq].clipped_after += words[0][
                    max(0,(far_bound - current_base)):]

                ## MODIFIED BEHAVIOUR:
                #  Take BETWEEN the rRNA and the tRNA.
                # This is "backwards" because the expected far bound is nearer when this happens.
                else:
                    seqs[clipping_seq].clipped_after += words[0][
                    max(0,(far_bound - current_base)):
                    max(0,near_bound - current_base)]

                ## This always happens
                # increment the counter
                current_base += len(words[0])

    for i in seqs:
        print(">"+i+"_cont_reg")
        seq = seqs[i].clipped_after+seqs[i].clipped_before
        if seqs[i].rev_comp:
            seq = reverseComp(seq)
        print(seq)

    seq_file.close()

# import-safety
if __name__ == '__main__':
    main()