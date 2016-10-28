class MafJunkContentQualityFilter:
    def __init__(self, acceptableNpercentage, junkCharacters):

        self.acceptableNpercentage = acceptableNpercentage
        self.junkchars = junkCharacters
        
    def accept(self, maf):
        junkchars = self.junkchars
        for i in range(maf.count()):                    
            ns = 0
            total = 0
            for j in maf.data(i):
                if j in junkchars:
                    ns += 1
                total += 1
                    
            if ns > self.acceptableNpercentage*total:
                return False
                
        return True
                    
                    
class MafIngroupJunkContentQualityFilter:
    def __init__(self, ingroup, acceptableNpercentage, junkCharacters):

        self.acceptableNpercentage = acceptableNpercentage
        self.junkchars = junkCharacters
        self.ingroup = ingroup
        
    def accept(self, maf):
        junkchars = self.junkchars
        for i in range(maf.count()):                    
            if maf.name(i) not in self.ingroup:
                continue
            ns = 0
            total = 0            
            for j in maf.data(i):
                if j in junkchars:
                    ns += 1
                total += 1
                    
            if ns > self.acceptableNpercentage*total:
                return False
                
        return True
                    
                    
