

class Annotation:
    def __init__(self, annotationFileName):
    
        annotationFile = open(annotationFileName)
        
        self.annotation = dict()
        
        for l in annotationFile:
	    lst = l.split()
	    self.annotation.setdefault(lst[0], []).append((int(lst[1]), int(lst[2])))
        	
        for key in self.annotation.iterkeys():
	    self.annotation[key].sort()
        
        annotationFile.close()
    		
    #binary search			
    def _contains(self, chrom, position, start, end):
    	if start >= end:
	    return False
    
    	m = start + (end - start)/2
    
    	(a,b) = self.annotation[chrom][m]
    
    	if a <= position and position < b:
	    return True
    	
    	if position < a:
	    return self._contains(self.annotation[chrom], position, start, m)
    	else:
	    return self._contains(self.annotation[chrom], position, m+1, end)
    
    def contains(self, chrom, position):
    	if chrom not in self.annotation:
	    return False
    	return self._contains(chrom, position, 0, len(self.annotation[chrom]))
    
# class Annotation:
# 	def __init__(self, annotationFileName):
# 
# 		annotationFile = open(annotationFileName)
# 		
# 		self.annotation = dict()
# 		
# 		for l in annotationFile:
# 			lst = l.split()
# 			self.annotation.setdefault(lst[0], []).append((int(lst[1]), int(lst[2])))
# 			
# 		for key in self.annotation.iterkeys():
# 			self.annotation[key].sort()
# 
# 		annotationFile.close()
# 			
# 	#binary search			
# 	def _contains(self, chrom, position, start, end):
# 		if start >= end:
# 			return False
# 
# 		m = start + (end - start)/2
# 	
# 		(a,b) = self.annotation[chrom][m]
# 	
# 		if a <= position and position < b:
# 			return True
# 		
# 		if position < a:
# 			return self._contains(self.annotation[chrom], position, start, m)
# 		else:
# 			return self._contains(self.annotation[chrom], position, m+1, end)
# 	
# 	def contains(self, chrom, position):
# 		if chrom not in self.annotation:
# 			return False
# 		return self._contains(chrom, position, 0, len(self.annotation[chrom]))
