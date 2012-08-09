#A class deciding whether 2 maf entries should merge
#this one takes ingroup into account
#contructed using
#IngroupMergeCriteria(
#parameter 1 : maximal distance acceptable for two maf entries to merge
#parameter 2 : list containing names of ingroups ie. ["hg18", bonobo"]
#)

#class deciding whether or not two mafs should be merged
class IngroupMergeCriteria:
	def __init__(self, maxGap, inGroup, maxSize=float('inf')):
		self.maxGap = maxGap
		self.inGroup = inGroup
		self.maxSize = maxSize
	
	#says true if they should be merged	
	def shouldMerge(self, maf1, maf2):

		if maf1.size(0) + maf2.size(0) > self.maxSize:
			return False
		
		for i in range(maf1.count()):
			assert maf1.name(i) == maf2.name(i)
			
			if maf1.name(i) not in self.inGroup:
				continue
				
			if maf1.strand(i) != maf2.strand(i):
				return False
			
			if maf1.chromosome(i) != maf2.chromosome(i):
				return False
				
			if maf2.start(i) - maf1.end(i) > self.maxGap:
				return False
		
		return True
