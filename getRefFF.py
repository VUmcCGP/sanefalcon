# Copyright (C) 2016 VU University Medical Center Amsterdam
# Author: Roy Straver (r.straver@vumc.nl)
#
# This file is part of SANEFALCON
# SANEFALCON is distributed under the following license:
# Attribution-NonCommercial-ShareAlike, CC BY-NC-SA (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
# This license is governed by Dutch law and this license is subject to the exclusive jurisdiction of the courts of the Netherlands.



import sys
import glob

def median(mylist):
    sorts = sorted(mylist)
    length = len(sorts)
    if not length % 2:
        return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
    return float(sorts[length / 2])
    
inFiles=glob.glob(sys.argv[1]+".*.start.*")

# Prepare bins per chromosome by their sizes
binSize = 50000
chrSizes = [16571,249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566]
chrArrays=[]
for c in chrSizes:
	chrArrays.append([0]*(c/binSize+1))
	#print c,len(chrArrays[-1])

#exit()

for inFile in inFiles:
	chrom = inFile.split(".")[-3]
	numChrom = 0
	# Turn X into the right chromosome number in our arrays
	if chrom == 'X':
		numChrom = 23
	else:
		numChrom = int(chrom)
	curChrArr = chrArrays[numChrom]
	with open(inFile) as reader:
		for line in reader:
			pos = int(line)/binSize
			curChrArr[pos] += 1
		
autosomals=[]
for chrom in chrArrays[1:23]:
	autosomals.extend(chrom)
medA=median(autosomals)
medX=median(chrArrays[23])

print sys.argv[1],sys.argv[1].split("/")[-1],((medA-medX)/medA*2)*100,medA,medX
		
'''
chrM	16571
chr1	249250621
chr2	243199373
chr3	198022430
chr4	191154276
chr5	180915260
chr6	171115067
chr7	159138663
chr8	146364022
chr9	141213431
chr10	135534747
chr11	135006516
chr12	133851895
chr13	115169878
chr14	107349540
chr15	102531392
chr16	90354753
chr17	81195210
chr18	78077248
chr19	59128983
chr20	63025520
chr21	48129895
chr22	51304566
chrX	155270560
chrY	59373566
'''

