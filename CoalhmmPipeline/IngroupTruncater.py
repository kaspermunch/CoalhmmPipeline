from MAF import *

class IngroupTruncater:
    def __init__(self, ingroup):
        self.ingroup = ingroup
        
    def truncate(self, maf):
        out = "a score=%f\n" % maf.score()
        
        lt = 0
        rt = 0
        for i  in range(len(maf.data(0))):
            ns = 0
            for j in range(maf.count()):
                if maf.name(j) not in self.ingroup:
                    continue
                    
                if maf.data(j)[i] in "N-":
                    ns = ns +1
                    
            if ns == len(self.ingroup):
                lt = lt +1
            else:
                break
                
        for i in range(len(maf.data(0))-1, 0, -1):
            ns = 0
            for j in range(maf.count()):
                assert len(maf.data(j)) == len(maf.data(0))
                if maf.name(j) not in self.ingroup:
                    continue
                    
                if maf.data(j)[i] in "N-":
                    ns = ns +1
                    
            if ns == len(self.ingroup):
                rt = rt +1
            else:
                break
                
        if (lt + rt) >= len(maf.data(0)):
            return None
            
        for j in range(maf.count()):
            assert len(maf.data(j)) == len(maf.data(0))
            
            
            start = maf.start(j)
            size = 0 #recalculated later
            strand = 1
            
            #trim data
            data = maf.data(j)[lt:len(maf.data(0))-rt]
            
            if maf.strand(j) == '-':
                strand = -1
                
            #adjust start
            for c in maf.data(j)[0:lt]:
                if c not in "-":
                    start = start + strand
                    
            #calculate size
            for c in data:
                if c not in "-":
                    size = size +1
                    
            out += "s %s.%s %i %i %s %i %s\n" % (maf.name(j), maf.chromosome(j), start, size, maf.strand(j), maf.srcLength(j), data)
            
        return MAF(out)
            
                    
