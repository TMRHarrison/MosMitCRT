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

###########################################################################
##     !!!     THIS ONLY WORKS IF IT'S IN THE FORMAT I EXPECT     !!!    ##
##     Which is: start -> Ileu tRNA -> other shit -> 12S rRNA -> end     ##
##                                                                       ##
##              If it isn't formatted like this, it'll grab              ##
##                  everything BUT the control sequence                  ##
###########################################################################


import sys
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "The file to be worked on")

args = parser.parse_args()

file = open(args.input)

bound1Pattern = re.compile("product=mtRNA-Ile\(...\)")
bound2Pattern = re.compile("product=12S ribosomal RNA")

currentSeq = ""

clipUntil = 0
clipAfter = 0

clipCounter = -1

clippedBefore = ""
clippedAfter = ""

for i in file:
	line = i[:-1]
	words = line.split(None, 8)

	# find the species name if you don't have one
	if currentSeq == "" and words[0] == "##sequence-region":
		currentSeq = words[1]
		print(currentSeq)

	# find the indices to start and stop at
	if words[0] == currentSeq:
		if bound1Pattern.search(words[8]):
			clipUntil = int(words[3]) - 1
		if bound2Pattern.search(words[8]):
			clipAfter = int(words[4])

	# stop when you reach a different sequence
	if words[0][:1] == ">":
		clipCounter = -1

	# clip when on the right sequence
	if clipCounter >= 0:
		# clip from the start until whatever
		clippedBefore += words[0][:clipUntil]
		clipUntil -= len(words[0])

		# Clip after getting to a certain point
		clippedAfter += words[0][max(0,(clipAfter - clipCounter)):]

		# increment the counter
		clipCounter += len(words[0])

	if words[0] == ">"+currentSeq:
		print("clipUntil: "+str(clipUntil))
		clipCounter = 0

	## old debug shit
	##for j in words:
	##	print(j)

print(currentSeq)
print(clippedAfter+"  "+clippedBefore)

file.close()