#!/usr/bin/python

# HERE'S THE PLAN
#
# Step 1: clip everything before the isoleucine into a new file
#
# Step 2: Clip everything after the small rRNA into said file
#
# I'll probably have to go back and figure out what to do if the rRNA
# appears before the Ileu but I'll figure that out later.
#
# New plan: Look at tRNA annotation. Look at rRNA annotation.
# Order as necessary:
# get sequence before Ileu
# get sequence after small rRNA
#
# TODO:
#   Make the sequences print to files instead of being held in memory?
#   Make the files more memory-safe and close if something happens during runtime
#

###########################################################################
##     !!!     THIS ONLY WORKS IF IT'S IN THE FORMAT I EXPECT     !!!    ##
##  Which is: start -> Ileu tRNA -> genes,rRNA,etc. -> 12S rRNA -> end   ##
##                                                                       ##
##              If it isn't formatted like this, it'll grab              ##
##                  everything BUT the control sequence                  ##
###########################################################################


import argparse
import re

def getParams():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help = "The file to be worked on")

    return parser.parse_args()

class Seq:
    name = ""
    nearBound1 = -1
    farBound1 = -1
    nearBound2 = -1
    farBound2 = -1
    clipped_before = ""
    clipped_after = ""
    innerFlag = False
    revComp = False

# check if the features has been flipped from the expected orientation (5'-B---A-3' instead of 5'-A---B-3')
def checkFlip(curSeq):
    ## safety for sequences when they're both found, and the 12S rRNA is before the ileu-tRNA
    if (curSeq.nearBound1 != -1
    and curSeq.farBound1 != -1
    and curSeq.farBound1 < curSeq.nearBound1):
        curSeq.innerFlag = (not curSeq.innerFlag) # take from between instead of the tails
        # it has to be a toggle instead of an assignment because of weirdnesss with rotated reverse complements.

# gets the reverse complement, don't worry about single letter variable names :)
def reverseComp(seq):
    d = {"A":"T", "T":"A", "G":"C", "C":"G"}
    o = ""
    for i in seq:
        o += d[i]
    return o[::-1]

def main():

    args = getParams()

    seq_file = open(args.input)

    bound1_pattern = re.compile("product=mtRNA-Ile\(...\)") ## escaped regex
    bound2_pattern = re.compile("product=12S ribosomal RNA")

    seqs = {}

    clip_counter = -1
    current_clip = None

    for i in seq_file:
        line = i[:-1] # clip the newline off
        words = line.split(None, 8) # split by whitespace

        # find the species name if you don't have one
        if words[0] == "##sequence-region":
            # Uses the sequence name as a reference so we don't have to loop
            # through the whole set every time we encounter something
            seqs[words[1]] = Seq()

        # find the indices to start and stop at
        if words[0] in seqs:
            curSeq = seqs[words[0]]


            ## only use the first tRNA-Ile
            if (curSeq.nearBound1 < 0
            and bound1_pattern.search(words[8])):
                curSeq.nearBound1 = int(words[3]) - 1
                curSeq.nearBound2 = int(words[4])
                checkFlip(curSeq)
            if bound2_pattern.search(words[8]):
                curSeq.farBound1 = int(words[3]) - 1
                curSeq.farBound2 = int(words[4])
                # operating under the assumption that the 12S rRNA is on the negative strand
                # so if it isn't, we assume that we're looking at the reverse complement of what we're normally looking at
                if words[6] == "+":
                    curSeq.revComp = True
                    # ok this part is kinda fucked
                    # we have to set it to be "flipped" because of the orientation of the sequence.
                    curSeq.innerFlag = True

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

                checkFlip(curSeq)

        # when you reach a different sequence, start processing it
        if words[0][:1] == ">":
            clip_counter = 0
            current_clip = words[0][1:]

        # if the line doesn't start with >, and clipCounter >= 0
        # then you start
        elif clip_counter >= 0:
            # if we're looking at the reverse complement, we need to flip the boundaries.
            if seqs[current_clip].revComp:
                nearBound = seqs[current_clip].farBound1
                farBound = seqs[current_clip].nearBound2
            else:
                nearBound = seqs[current_clip].nearBound1
                farBound = seqs[current_clip].farBound2



            ## DEFAULT BEHAVIOUR:
            #  take from before the tRNA and then after the rRNA
            # clip from the start until whenever the tRNA begins
            if not seqs[current_clip].innerFlag:
                seqs[current_clip].clipped_before += words[0][
                :max(0,nearBound - clip_counter)]

                # Clip after getting to the end of the rRNA
                seqs[current_clip].clipped_after += words[0][
                max(0,(farBound - clip_counter)):]

            ## MODIFIED BEHAVIOUR:
            #  Take BETWEEN the rRNA and the tRNA.
            # This is "backwards" because the expected far bound is nearer when this happens.
            else:
                seqs[current_clip].clipped_after += words[0][
                max(0,(farBound - clip_counter)):
                max(0,nearBound - clip_counter)]

            ## This always happens
            # increment the counter
            clip_counter += len(words[0])

    for i in seqs:
        print(">"+i+"_cont_reg")
        seq = seqs[i].clipped_after+seqs[i].clipped_before
        if seqs[i].revComp:
            seq = reverseComp(seq)
        print(seq)

    seq_file.close()

# import-safety
if __name__ == '__main__':
    main()