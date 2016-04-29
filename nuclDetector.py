# Copyright (C) 2016 VU University Medical Center Amsterdam
# Author: Roy Straver (r.straver@vumc.nl)
#
# This file is part of SANEFALCON
# SANEFALCON is distributed under the following license:
# Attribution-NonCommercial-ShareAlike, CC BY-NC-SA (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
# This license is governed by Dutch law and this license is subject to the exclusive jurisdiction of the courts of the Netherlands.



import sys
import operator

areaSize=73 #60
sideSize=20 #25
padding=areaSize+sideSize
outerLen=2*sideSize#len(left)+len(right)
innerLen=areaSize*2+1#len(window)


if len(sys.argv) > 3:
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt

def plotRegion(start,end,area,bpScores,centers):
	if len(centers) < 6 or len(area) < 2000:
		return
	
	plt.clf()
	
	left=len(area)/2-1000
	right=len(area)/2+1000
	
	smooth=[]
	for i in range(len(area)):
		smooth.append(sum(area[max(0,i-10):i+10])/20.)
	print smooth
	print centers
	for center in centers:
		#plt.axvline(center[0]-start-74, color="gray")
		#plt.axvline(center[0]-start+74, color="gray")
		plt.axvspan(center[0]-73, center[0]+73, color='gray', alpha=0.3)
		plt.axvline(center[0], color="black")
	plt.plot(smooth)
	plt.plot(bpScores,"red")
	plt.title("Smoothed Read Start Count (chr21:"+str(start+left)+"-"+str(start+right)+")")
	plt.xlabel("Relative BP position")
	plt.ylabel("Count or Score")
	plt.xlim([left,right])
	plt.ylim([0,10])
	fig = matplotlib.pyplot.gcf()
	fig.set_size_inches(16,2)
	plt.savefig(sys.argv[3]+"_"+str(start)+"-"+str(end)+'.plot.pdf', dpi=100)


def flush(area,endPoint):
	# Sliding window settings, area is half the nucleosome size
	bins=[0]*padding+area+[0]*padding
	areaLen=len(area)
	startPoint=endPoint-areaLen+1 # TODO check +1?

	# Let's move out
	bpScores=[]
	extra=[]
	peaks=[]

	# Slide  a window of 20\147/20 over the region and score positions
	for i in range(padding,areaLen+padding):
		window=bins[i-areaSize:i+areaSize+1]
		
		leftVal=sum(bins[i-areaSize-sideSize:i-areaSize])
		rightVal=sum(bins[i+areaSize+1:i+areaSize+sideSize+1])
		innerVal=sum(window)
		outerVal=leftVal+rightVal#min(leftVal,rightVal)
		#if leftVal > 25 and rightVal > 25:
			# Solve zero devision and score this pos twice as high
		#	innerVal=max(innerVal,0.5)
		score = ((outerVal+1)/float(outerLen))/((innerVal+1)/float(innerLen))
		bpScores.append(score)
		extra.append([leftVal,innerVal,rightVal])
		#else:
		#	bpScores.append(0)
	
	# Pick best scoring positions in the region until no positive values exist
	def findCenters(start,end):
		if end-start < 1:
			return []
		maxIndex, bpMax = max(enumerate(bpScores[start:end]), key=operator.itemgetter(1))
		if bpMax > 1:
		
			# If the max is part of a flat line take it's center position instead
			tmpIndex=maxIndex
			while tmpIndex < ( end-start-1 ) and bpScores[ tmpIndex + 1 ] == bpMax:
				tmpIndex += 1
			maxIndex  = ( maxIndex + tmpIndex ) / 2
			
			left      = start + maxIndex - innerLen
			right     = start + maxIndex + innerLen + 1
			leftList  = findCenters( start , left )
			rightList = findCenters( right , end  )
			thisList  = [ start + maxIndex , bpMax ]
			return leftList + [ thisList ] + rightList
		return []
		
	newCenters=findCenters(0,areaLen)
	if newCenters!=[] and newCenters[-1]!=[] and newCenters[-1][0] > len(extra):
		print newCenters[-1],len(extra),len(bpScores),bpScores[-5:]
	allNucl.extend([[x[0]+startPoint]+x[1:]+extra[x[0]] for x in newCenters])
	
	
	
	if len(sys.argv) > 3:
		plotRegion(startPoint,endPoint,area,bpScores,newCenters)
	
curArea=[0]
lastPos=0
maxDist=190 # A little over our sliding window size
inFile = open(sys.argv[1], 'rb')
allNucl=[]

for line in inFile:
	splitLine=line.split()
	position=int(splitLine[0])
	distance = position - lastPos
	
	#if position > 17500000:
	#	break

	# So far away no nucleosome can be detected in between, finish this region 
	if distance > maxDist:
		flush(curArea,lastPos)
		curArea=[1]
	# Add read to current region
	else:
		curArea+=[0 for x in range(distance)]
		#curArea.extend([0 for x in range(distance+1000)])
		curArea[-1]+=1
	lastPos = position
flush(curArea,lastPos)
inFile.close()

# Dump results to a file
with open(sys.argv[2], "wb") as output_file:
	for nucl in allNucl:
		#output_file.write(str(nucl[0]) + "\t" + str(nucl[1]) + "\n")
		output_file.write("\t".join([str(x) for x in nucl])+"\n")
		
		
'''
if distance > maxDist:
	
	if lastPos >= int(sys.argv[4]):
		print "End of this process reached:",sys.argv[4],lastPos
		inFile.close()
		quit()
	elif lastPos >= int(sys.argv[3]):
		flush(curArea,lastPos)
	else:
		print "Start of this process not yet reached:",sys.argv[3],lastPos
	curArea=[1]
else:
'''
