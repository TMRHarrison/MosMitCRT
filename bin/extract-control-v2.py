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
    clip_until = -1
    clip_after = -1
    clipped_before = ""
    clipped_after = ""

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
            ## only use the first tRNA-Ile
            if (seqs[words[0]].clip_until < 0
            and bound1_pattern.search(words[8])):
                seqs[words[0]].clip_until = int(words[3]) - 1
            if bound2_pattern.search(words[8]):
                seqs[words[0]].clip_after = int(words[4]) + 1

        # when you reach a different sequence, start processing it
        if words[0][:1] == ">":
            clip_counter = 0
            current_clip = words[0][1:]

        # if the line doesn't start with >, and clipCounter >= 0
        # then you start
        elif clip_counter >= 0:

            # clip from the start until whenever the tRNA begins
            seqs[current_clip].clipped_before += words[0][
            :max(0,seqs[current_clip].clip_until - clip_counter)]

            # Clip after getting to the end of the rRNA
            seqs[current_clip].clipped_after += words[0][
            max(0,(seqs[current_clip].clip_after - clip_counter)):]

            # increment the counter
            clip_counter += len(words[0])

    for i in seqs:
        print(">"+i+"_cont_reg")
        print(seqs[i].clipped_after+seqs[i].clipped_before)

    seq_file.close()

# import-safety
if __name__ == '__main__':
    main()