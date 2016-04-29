# Copyright (C) 2016 VU University Medical Center Amsterdam
# Author: Roy Straver (r.straver@vumc.nl)
#
# This file is part of SANEFALCON
# SANEFALCON is distributed under the following license:
# Attribution-NonCommercial-ShareAlike, CC BY-NC-SA (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
# This license is governed by Dutch law and this license is subject to the exclusive jurisdiction of the courts of the Netherlands.



# 1 Peak file as reference for nucleosome positions
# 2 Read start pos file for a sample
# 3 Reverse: 0 or 1
# 4 Output basename

import sys
import math

minSide=25
minCenter=100

areaSize=73. #60
sideSize=20. #25

filtOut=0

def nuclFilt(prop):
	global filtOut
	#if min(prop[0],prop[2]) > minSide and prop[1] > minCenter and min(prop[0],prop[2])/sideSize > prop[1]/147.: # ex8
	#if min(prop[0],prop[2])/sideSize > 1.5*prop[1]/147. and prop[1] > minCenter: #ex9
	if min(prop[0],prop[2])/sideSize > prop[1]/147.: # and prop[1] > minCenter: #ex9
		return True
	filtOut+=1
	return False

def loadNucl(nuclLine):
	#return [int(nuclLine[0]),(min(float(nuclLine[2]),float(nuclLine[4]))/sideSize)/(float(nuclLine[3])/147)]
	return [int(nuclLine[0]),float(nuclLine[1])]
	
#peaks = [[int(line.split()[0]),float(line.split()[1])] for line in open(sys.argv[1]) if nuclFilt([int(x) for x in line.split()[2:]])]
peaks = [loadNucl(line.split()) for line in open(sys.argv[1]) if nuclFilt([int(x) for x in line.split()[2:]])]
#peaks = [[int(line.split()[0]),float(line.split()[1])] for line in open(sys.argv[1])]
#peaks = [[int(line.split()[0]),1] for line in open(sys.argv[1]) if float(line.split()[1])>1.5]
reads = [int(line) for line in open(sys.argv[2])]
shift = int(sys.argv[3])
print len(peaks),filtOut

#for nuc in peaks[:10]:
#	print nuc

#exit()

rev=False
if shift != 0:
	rev=True
	peaks.reverse()
	reads.reverse()
	# firx for 100bp pe data
	#reads=[x+50 for x in reads]
	#peaks=[[x[0]-50,x[1]] for x in peaks]
#else:
	#reads=[x+50 for x in reads]
	#peaks=[[x[0]+50,x[1]] for x in peaks]
peaks.append([-1,0])

maxDist=147
sumPeak=[0. for x in range(maxDist)]
read = reads[0]#+shift
j=0
nuclHit=[]

#print peaks[:10],reads[:10]

if rev:
	for i,peakPair in enumerate(peaks):
		peak=peakPair[0]
		peakWeight=peakPair[1]
		thisPeak=[0. for x in range(maxDist)]
		while read >= peak:
			if read < peak+maxDist:
				thisPeak[read-peak] += 1
			j+=1
			if j >= len(reads):
				break
			read=reads[j]#+shift
		thisSum=float(sum(thisPeak))
		if thisSum == 0:
			continue
			
		nuclHit.append(peak)
		#thisPeakNorm=[x/thisSum for x in thisPeak]
		#thisPeak=[x*peakWeight for x in thisPeak]
		sumPeak=[sumPeak[x]+thisPeak[x] for x in range(len(thisPeak))]
	nuclHit.reverse()

else:
	for i,peakPair in enumerate(peaks):
		peak=peakPair[0]
		peakWeight=peakPair[1]
		thisPeak=[0. for x in range(maxDist)]
		while read <= peak:
			if read > peak-maxDist:
				thisPeak[peak-read] += 1
			j+=1
			if j >= len(reads):
				break
			read=reads[j]#+shift

		thisSum=float(sum(thisPeak))
		if thisSum == 0:
			continue
		
		nuclHit.append(peak)
		#thisPeakNorm=[x/thisSum for x in thisPeak]
		#thisPeak=[x*peakWeight for x in thisPeak]
		sumPeak=[sumPeak[x]+thisPeak[x] for x in range(len(thisPeak))]

#print nuclHit
	
with open(sys.argv[4]+".np",'w') as outPeaks:
	outPeaks.write(",".join([str(x) for x in sumPeak]))
exit()
with open(sys.argv[4]+".nucl",'w') as outNucl:
	for nuclUsed in nuclHit:
		outNucl.write(str(nuclUsed)+'\n')

