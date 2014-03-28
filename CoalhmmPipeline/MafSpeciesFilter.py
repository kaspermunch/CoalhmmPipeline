
class MafSpeciesFilter:
    def __init__(self, species):
        self.species = species

    def inplace(self, maf):
        """
        Filter MAF in place leaving only specified species
        """
        for idx in reversed(range(maf.count())):
            if maf.name(idx) not in self.species:
                del maf.alignment[idx]        
            
