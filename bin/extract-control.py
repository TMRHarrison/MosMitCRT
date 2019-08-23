#!/usr/bin/env python

#   Make the sequences print to files instead of being held in memory?
#   Extend the system to clip any control region, given command line arguments for the annotation boundaries

## Naming:
# variable_names
# functionNames
# Classes

# command line arguments, regular expressions
import argparse
import re

def getParams():
    """
    Returns the command line arguments.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="The file to be worked on. GFF format with a ##FASTA section.")

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

    def checkFlip(self):
        """
        toggles inner_flag if the features has been flipped from the expected orientation (5'-B---A-3' instead of 5'-A---B-3').
        """

        ## safety for sequences when they're both found, and the END sequence is after
        if (self.start_anno_5_prime != -1
                and self.end_anno_5_prime != -1
                and self.end_anno_5_prime < self.start_anno_5_prime):
            self.inner_flag = (not self.inner_flag) # take from between instead of the tails
            # it has to be a toggle instead of an assignment because of weirdnesss with rotated reverse complements.

    def getBounds(self, words, bound_start_pattern, bound_end_pattern):
        """
        Grabs all the nt in the control sequence and puts them in clipped_before an clipped_after.
        """

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
        if (self.start_anno_5_prime < 0
                and bound_start_pattern.search(words[8])):
            self.start_anno_5_prime = int(words[3]) - 1
            self.start_anno_3_prime = int(words[4])
            self.checkFlip()
        if bound_end_pattern.search(words[8]):
            self.end_anno_5_prime = int(words[3]) - 1
            self.end_anno_3_prime = int(words[4])
            # operating under the assumption that the 12S rRNA is on the negative strand
            # so if it isn't, we assume that we're looking at the reverse complement of what we're normally looking at
            if words[6] == "+":
                self.rev_comp = True
                # we have to set inner flag to be "flipped" because of the orientation of the sequence.
                self.inner_flag = True

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
            self.checkFlip()

    def parseSeq(self, words, current_base):
        """
        Clips the relevant sequence data from sequences
        """

        # if the line doesn't start with >, and clipCounter >= 0
        # then you start
        if (current_base >= 0) and (words[0][:1] != ">"):
            # if we're looking at the reverse complement, we need to flip the boundaries.
            if self.rev_comp:
                near_bound = self.end_anno_5_prime
                far_bound = self.start_anno_3_prime
            else:
                near_bound = self.start_anno_5_prime
                far_bound = self.end_anno_3_prime

            ## DEFAULT BEHAVIOUR:
            #  take from before the tRNA and then after the rRNA
            # clip from the start until whenever the tRNA begins
            if not self.inner_flag:
                self.clipped_before += words[0][
                    :max(0, near_bound - current_base)]

                # Clip after getting to the end of the rRNA
                self.clipped_after += words[0][
                    max(0, (far_bound - current_base)):]

            ## MODIFIED BEHAVIOUR:
            #  Take BETWEEN the rRNA and the tRNA.
            # This is "backwards" because the expected far bound is nearer when this happens.
            else:
                self.clipped_after += words[0][
                    max(0, (far_bound - current_base)):
                    max(0, near_bound - current_base)]

            ## This always happens
            # return the number of nt parsed so we can incrment the counter
            return len(words[0])
        # if we find nothing, just return 0.
        return 0

# biopython can also do this, but I feel like it's easier to do it this way if I just need the one thing.
# makes a translation dictionary in case of reverse compliments
NT_ALL  = "AGCTURYKMBVDH" # everything else will get ignored -> '-', 'N', "S", "W" get "complemented" to themselves
NT_COMP = "TCGAAYRMKVBHD"
# build the dictionary and store it
NT_COMP_DICT = str.maketrans(NT_ALL+NT_ALL.lower(), NT_COMP+NT_COMP.lower())

# we can just pop these out since we don't need them any more and they're just polluting the namespace at this point
del NT_ALL, NT_COMP

# Uses the library defined just above to get the other strand, then reverses the string
def reverseComp(seq):
    """
    Gives the reverse complement of a sequence input.
    """
    return seq.translate(NT_COMP_DICT)[::-1]

def main():
    """
    Does the whole thing
    I don't know what you want me to say
    """
    args = getParams()

    # dictionary of sequences indexed by name
    seqs = {}

    with open(args.input, 'r') as seq_file:

        # The annotation you want the cut to START at
        bound_start_pattern = re.compile(r'product=mtRNA-Ile\(...\)') ## escaped regex
        # The annotation the cut ENDS at
        bound_end_pattern = re.compile(r'product=12S ribosomal RNA')

        # current base position being parsed
        current_base = -1

        # name of the sequence being clipped
        clipping_seq = None

        for line in seq_file:
            line = line[:-1] # clip the newline off
            words = line.split(None, 8) # split by whitespace

            if words[0] == "##sequence-region":
                # adds a new sequence to the list of sequences when it finds a sequence-region tag.
                seqs[words[1]] = SeqInfo()
                # Uses the sequence name as a reference so we don't have to loop
                # through the whole set every time we encounter something

            if words[0] in seqs:
                seqs[words[0]].getBounds(words, bound_start_pattern, bound_end_pattern)

            if words[0][:1] == ">":
                current_base = 0
                clipping_seq = words[0][1:]

            if clipping_seq in seqs:
                current_base += seqs[clipping_seq].parseSeq(words, current_base)

    for i in seqs:
        print(">"+i+"_cont_reg")
        seq = seqs[i].clipped_after+seqs[i].clipped_before
        if seqs[i].rev_comp:
            seq = reverseComp(seq)
        print(seq)

# import-safety
if __name__ == '__main__':
    main()