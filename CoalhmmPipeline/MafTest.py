

class MafTest(object):

    def __init__(self, speciesList):
        self.speciesList = sorted(speciesList)

    def test(self, maf):
        return self.speciesList == sorted([maf.name(i) for i in range(maf.count())])

        

        
