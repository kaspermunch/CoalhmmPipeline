class MafNContentQualityFilter:
    def __init__(self, acceptableNpercentage):

        self.acceptableNpercentage = acceptableNpercentage
        
    def accept(self, maf):
        for i in range(maf.count()):                    
            ns = 0
            total = 0
            for j in maf.data(i):
                if j == 'N':
                    ns += 1
                    total += 1
                elif j in "ACGTacgt":
                    total = total +1
                    
            if ns > self.acceptableNpercentage*total:
                return False
                
        return True
                    
                    
