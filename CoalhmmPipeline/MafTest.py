

class MafTest(object):

    def __init__(self, speciesList):
        self.speciesList = speciesList

    def test(self, maf):
        return self.speciesList == [maf.name(i) for i in range(maf.count())]

        

        
