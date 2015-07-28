

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
    
    	m = start + (end - start)/float(2)
    
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

"""
    def contains(self, chrom, position):
    	if chrom not in self.annotation:
	    return False
    	return self._contains(self.annotation[chrom], position, startIdx=1, endIdx=2)

    def _contains(lst, pos, startIdx, endIdx, lo=0, hi=None):
        if hi is None:
            hi = len(lst)
        while lo < hi:
            mid = (lo+hi)//2
            start, end = lst[mid][startIdx, endIdx] # <- change to extraction of pos value here
            if end <= pos:
                lo = mid+1
            elif start > pos: 
                hi = mid
            else:
                return True
        return False
"""

#     def binary_search(a, x, lo=0, hi=None):
#         if hi is None:
#             hi = len(a)
#         while lo < hi:
#             mid = (lo+hi)//2
#             midval = a[mid] # <- change to extraction of pos value here
#             if midval < x:
#                 lo = mid+1
#             elif midval > x: 
#                 hi = mid
#             else:
#                 return mid
#         return -1
#     # depending on whether the upper bound is inclusive or exclusive.
#     # you can change hi = mid to hi = mid-1 and hi = len(a) to hi = len(a)-1 and while lo < hi: to while lo <= hi,
#     # and it would be equivalently correct 



    
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
