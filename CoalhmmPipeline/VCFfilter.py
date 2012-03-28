
class VCFfilter():

    def __init__(self, filters):    
        self.filters = [eval(f) for f in filters]

    def accept(self, snp):
        return all(f(snp) for f in self.filters)

    def no(self, snp):
        return True

    def Q30(self, snp):
        """
        Filter on phred score 30
        """
        pass

    def Q50(self, snp):
        """
        Filter on phred score 50
        """
        pass

    def C20(self, snp):
        """
        Filter on coverage 20
        """
        pass
