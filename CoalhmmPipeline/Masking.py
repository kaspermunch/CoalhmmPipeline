import re, sys

from MAF import *
from array import array

class Masking(object):

    def __init__(self, annotation):
    	
        self.annotation = annotation
    
    def maskMafFile(self, referenceNameRegex, inputFileName, outputFileName):
    
    	inputFile = open(inputFileName)
        output = ''
        
    	referenceNameRegex = re.compile(referenceNameRegex)
    
    	for maf in MAFIterator(inputFile):
    	    #find index of the reference
    	    referenceIndex = None
    	    for i in range(len(maf.alignment)):
	        if referenceNameRegex.match(maf.alignment[i][1]):
		    referenceIndex = i
		    break
    	    assert referenceIndex is not None
    	    
    	    #find chromosome of the reference
    	    chrom = maf.alignment[referenceIndex][1].split(".")[1]

    	    mafStart = int(maf.alignment[referenceIndex][2])
    	    mafEnd = int(maf.alignment[referenceIndex][2]) + int(maf.alignment[referenceIndex][2])
    	    
    	    alignmentPosition = mafStart
    	    #convert strings into arrays
    	    for j in range(len(maf.alignment)):
	        maf.alignment[j][6] = array("c", maf.alignment[j][6])
    	    
    	    for i in xrange(len(maf.alignment[referenceIndex][6])):
	        currentChar = maf.alignment[referenceIndex][6][i]
    	    
	        if self.annotation.contains(chrom, alignmentPosition):
		    for j in xrange(len(maf.alignment)):
		        maf.alignment[j][6][i] = "N"
    	    
    	    	if currentChar != "-":
		    alignmentPosition += 1
    	    
    	    for j in range(len(maf.alignment)):
	        maf.alignment[j][6] = maf.alignment[j][6].tostring()
    	    
            output += str(maf)

    	outputFile = open(outputFileName, 'w')
        print >>outputFile, output

 
# class Masking(object):
# 
# 	def __init__(self, annotation):
# 		
# 		self.annotation = annotation
# 
# 	def maskMafFile(self, referenceNameRegex, inputFileName, outputFileName):
# 
# 		inputFile = open(inputFileName)
# 		outputFile = open(outputFileName, 'w')
# 
# 		referenceNameRegex = re.compile(referenceNameRegex)
# 
# 		for maf in MAFIterator(inputFile):
# 			#find index of the reference
# 			referenceIndex = None
# 			for i in range(len(maf.alignment)):
# 				if referenceNameRegex.match(maf.alignment[i][1]):
# 					referenceIndex = i
# 					break
# 			assert referenceIndex is not None
# 
# 			#find chromosome of the reference
# 			chrom = maf.alignment[referenceIndex][1].split(".")[1]
# 
# 			mafStart = int(maf.alignment[referenceIndex][2])
# 			mafEnd = int(maf.alignment[referenceIndex][2]) + int(maf.alignment[referenceIndex][2])
# 
# 			alignmentPosition = mafStart
# 			#convert strings into arrays
# 			for j in range(len(maf.alignment)):
# 				maf.alignment[j][6] = array("c", maf.alignment[j][6])
# 
# 			for i in xrange(len(maf.alignment[referenceIndex][6])):
# 				currentChar = maf.alignment[referenceIndex][6][i]
# 
# 				if self.annotation.contains(chrom, alignmentPosition):
# 					for j in xrange(len(maf.alignment)):
# 						maf.alignment[j][6][i] = "N"
# 
# 				if currentChar != "-":
# 					alignmentPosition += 1
# 
# 			for j in range(len(maf.alignment)):
# 				maf.alignment[j][6] = maf.alignment[j][6].tostring()
# 
# 			print >>outputFile, maf
