class MafGroupLengthQualityFilter:
    def __init__(self, minLengthIngroup, minLenghtOutgroup, ingroup):
        self.minLengthIngroup = minLengthIngroup
        self.minLengthOutgroup = minLenghtOutgroup
        self.ingroup = ingroup
        
    def accept(self, maf):
        for i in range(maf.count()):
            if maf.name(i) in self.ingroup:
                if maf.length(i) < self.minLengthIngroup:
                    return False
            else:
                if maf.length(i) < self.minLengthOutgroup:
                    return False                
        return True
                    
                    
            
        
