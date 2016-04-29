# Copyright (C) 2016 VU University Medical Center Amsterdam
# Author: Roy Straver (r.straver@vumc.nl)
#
# This file is part of SANEFALCON
# SANEFALCON is distributed under the following license:
# Attribution-NonCommercial-ShareAlike, CC BY-NC-SA (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
# This license is governed by Dutch law and this license is subject to the exclusive jurisdiction of the courts of the Netherlands.



import sys
import pickle
import argparse

parser = argparse.ArgumentParser(description='Convert any stream of reads to a pickle file for WISECONDOR, defaults are set for the SAM format',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-binsize', type=int, default=1000000,
					help='binsize used for samples')
parser.add_argument('-retdist', type=int, default=4,
					help='maximum amount of base pairs difference between sequential reads to consider them part of the same tower')
parser.add_argument('-retthres', type=int, default=4,
					help='threshold for when a group of reads is considered a tower and will be removed')
parser.add_argument('-colchr', type=int, default=2,
					help='column containing chromosome, default is for sam format')
parser.add_argument('-colpos', type=int, default=3,
					help='column containing read start position, default is for sam format')

args = parser.parse_args()

binsize		= args.binsize
minShift 	= args.retdist
threshold	= args.retthres
chrColumn	= args.colchr
startColumn	= args.colpos

# Prepare the list of chromosomes
chromosomes = dict()
for chromosome in range(1,23):
	chromosomes[str(chromosome)] = [0]
chromosomes['X'] = [0]
chromosomes['Y'] = [0]


# Flush the current stack of reads
def flush(readBuff):
	global chromosomes
	stairSize = len(readBuff)
	if stairSize <= threshold or threshold < 0:	
		for line in fullBuff:
			print line.strip()
				
						
prevWords = ['0'] * 10
readBuff = []
fullBuff = []
for line in sys.stdin:
	curWords = line.split()

	# Not ndup, flush and start new stair
	if not((curWords[chrColumn] == prevWords[chrColumn]) and (minShift >= (int(curWords[startColumn])-int(prevWords[startColumn])))):
		flush(readBuff)
		readBuff = []
		fullBuff = []

	# Ignore reads starting at the exact same position: likely PCR dups
	if len(readBuff) ==0 or int(curWords[startColumn]) != readBuff[-1][1]:
		# Normal ndups will be appended here
		readBuff.append([curWords[chrColumn],int(curWords[startColumn])])
		fullBuff.append(line)
	
	prevWords = curWords
	prevLine = line

# Flush after we're done
flush(readBuff)
