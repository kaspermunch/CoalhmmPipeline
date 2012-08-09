#generic chunkenizer using 3 objects for splitting, criteria and output		
#Construct using
#Chunkenizer(
#parameter 1 : splitter - object potetially manipulating maf entries shortening or even splitting
#parameter 2 : merge criteria - a decider choosen whether or not two maf entries should be merged
#parameter 3 : output - an object responsible outputting the mafs in the chunks described by the criteria
#)
#to execute just call
#chunkenizerObject.chunkenize(
#parameter 1 : string containing the filename for the maf that should be chunkenized
#)
#WARNING: Use instances only once! depending parameter objects data may be corrupted
#TODO: consider making a use once or just a function
from MAF import *
from collections import deque

class Chunkenizer:
	def __init__(self, splitter, mergeCriteria, mafTest, mafQualityFilters, chunkQualityFilters, truncater):
		self.mafTest = mafTest
		self.splitter = splitter
		self.mergeCriteria = mergeCriteria
		self.mafQualityFilters = mafQualityFilters
		self.chunkQualityFilters = chunkQualityFilters
		self.truncater = truncater
	
	def acceptMaf(self, maf):
		for mqf in self.mafQualityFilters:
			if not mqf.accept(maf): 
				return False
		
		return True
		
	def acceptChunk(self, mafs):
		for cqf in self.chunkQualityFilters:
			if not cqf.accept(mafs): 
				return False
		
		return True
		
	def chunkenize(self, inputMAF, output):

		toMerge = [] #list of mafs to merge
		q = deque() #queue containing mafs which still need processing
		mafFile = open(inputMAF) 

# 		TEST_print = False
		
		for maf in MAFIterator(mafFile):

# 			print "INPUT"
# 			print maf
# 			assert all(maf.srcLength(i) > 0 for i in range(maf.count()))
# 			for i in range(len(maf.alignment)):
# 				_, TEST_name, TEST_start, TEST_length, _, _, _ = maf.alignment[i]
# 				if TEST_name == "Hsap" and TEST_start > 0.97e8:
# 					TEST_print = True
# 			if TEST_print:
# 				raw_input("press enter")
# 				print "Input MAF"
# 				print maf


			if not self.mafTest.test(maf):
#  				print maf
				continue

			#split maf and put it on the processing que
			for splitmaf in self.splitter.split(maf):
				#truncate the maf
				#print "before truncating\n", splitmaf, "\nafter"
				splitmaf = self.truncater.truncate(splitmaf)
				#print splitmaf
				#raw_input("press enter")    
				
				if splitmaf != None and self.acceptMaf(splitmaf):
					q.append(splitmaf)
			

			if len(q) == 0: #queue is empty get the next
			    continue

			"""if len(self.splitter.split(maf)) > 1:
				print "this\n", maf, "\nsplitinto"
				for splitmaf in self.splitter.split(maf):
					print splitmaf
					
				print "end split"
				raw_input("press enter")"""
				
			
			#if the tomerge list is empty then just add one maf from the queue
			#later code assumes toMerge is not empty
			if len(toMerge) == 0:
				toMerge.append(q.popleft())
			
			#now just flush the queue
			#for every element in the queue: test with the mergeCriteria against the last added maf in tomerge list
			#if they should merge: just add the element to the tomerge list
			#otherwise writeout the tomerge list, and then add.
			
			while len(q) > 0:
				e = q.popleft()
				last = toMerge[len(toMerge) -1]
				
				if self.mergeCriteria.shouldMerge(last, e):
					toMerge.append(e)
				else:
				    if self.acceptChunk(toMerge):
    					output.writeChunk(toMerge)

# 					for m in toMerge:
# 						print "OUTPUT"
# 						print m
# 						assert all(m.srcLength(i) > 0 for i in range(m.count()))
# 
# 					if TEST_print:
# 						for m in toMerge:
# 							print "Output MAF"
# 							print m
# 
					toMerge = [e]


		
		#since the tomerge list is never emty we need to write out the last chunk
		if len(toMerge) > 0 and self.acceptChunk(toMerge):
			output.writeChunk(toMerge)			
		
		#fullfill the output objects last wishes... if any
		output.finalize()



#example use		
#ingroup = ["hg18", "pantro2", "bonobo"]
#chunkenizer = Chunkenizer(Splitter(), MergeCriteria(kmerge, ingroup), Output("/tmp/chunks", ingroup, "/tmp/hdf.h5"))
#chunkenizer.chunkenize(inputFile)
