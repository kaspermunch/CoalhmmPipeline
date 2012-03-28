
from string import maketrans

from BufferedIterator import BufferedIterator

class VCF2MAF(object):

    def __init__(self, mafIterator, genomeIdx, refname, snpSpec):

        self.mafIterator = mafIterator
        self.genomeIdx = genomeIdx
        self.refname = refname
        self.snpIterators = dict()
        for species in snpSpec:
            self.snpIterators[species] = BufferedIterator(snpSpec[species])
        self.speciesNames = sorted(snpSpec.keys())

    def self._getSubstitutions(species, chrom, start, end):
        snpList = list()
        while True:
            snp = self.snpIterators[species].next()
            if snp.pos < start:
                continue
            elif snp.pos >= end:
                self.snpIterators[species].putback(snp)
                break
            else:
                snpList.append(snp)
        substitutions = [(snp.pos - start, snp.data["ALT"], snp.data["REF"]) for snp in snpList]
        return substitutions

    def writeMaf(self, outputFileName):

        outputFile = open(outputFile)

        complement = maketrans("AGTC", "TCAG")
        
        for maf in self.mafIterator:

            refIdx = None
            for i in range(len(maf.alignment)):
                    if referenceNameRegex.match(maf.alignment[i][1]):
                            refIdx = i
                            break
            assert refIdx is not None

            chrom = maf.chromosome(refIdx)
            srcLength, start, length = maf.srclength(refIdx), maf.start(refIdx), maf.length(refIdx)

            if strand == '-' :
                start = srcLength - start
            
            seq = self.genomeIdx.get(name, start, length)          

            if strand == '-':
                seq = seq.translate(complement)[::-1]

            speciesSeq = seq
            for species in self.speciesNames:
                for idx, subst, ref in self._getSubstitutions(species, chrom, start, end):
                    assert speciesSeq[idx] == ref
                    speciesSeq[idx] = subst

            maf.setSeq(idx, speciesSeq)

            # tweek annotation to reflect that this is not a real aligmne

            print >>outputFile, newMaf


    outputFile.close()
