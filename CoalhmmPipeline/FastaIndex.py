
import sys, os, string
import cPickle as pickle

class ProgressBar(object):

    def __init__(self, minValue = 0, maxValue = 100, totalWidth=60):
        self.progBar = "[]"   # This holds the progress bar string
        self.min = minValue
        self.max = maxValue
        self.span = maxValue - minValue
        self.width = totalWidth
        self.amount = 0       # When amount == max, we are 100% done 
        self.updateAmount(0)  # Build progress bar string

    def updateAmount(self, newAmount = 0):
        """ Update the progress bar with the new amount (with min and max
            values set at initialization; if it is over or under, it takes the
            min or max value as a default. """
        if newAmount < self.min: newAmount = self.min
        if newAmount > self.max: newAmount = self.max
        self.amount = newAmount

        # Figure out the new percent done, round to an integer
        diffFromMin = float(self.amount - self.min)
        percentDone = (diffFromMin / float(self.span)) * 100.0
        percentDone = int(round(percentDone))

        # Figure out how many hash bars the percentage should be
        allFull = self.width - 2
        numHashes = (percentDone / 100.0) * allFull
        numHashes = int(round(numHashes))

        # Build a progress bar with an arrow of equal signs; special cases for
        # empty and full
        if numHashes == 0:
            self.progBar = "[>%s]" % (' '*(allFull-1))
        elif numHashes == allFull:
            self.progBar = "[%s]" % ('='*allFull)
        else:
            self.progBar = "[%s>%s]" % ('='*(numHashes-1),
                                        ' '*(allFull-numHashes))

        # figure out where to put the percentage, roughly centered
        percentPlace = (len(self.progBar) / 2) - len(str(percentDone)) 
        percentString = str(percentDone) + "%"

        # slice the percentage into the bar
        self.progBar = ''.join([self.progBar[0:percentPlace], percentString,
                                self.progBar[percentPlace+len(percentString):]
                                ])

    def __str__(self):
        return str(self.progBar)


    def progressHook(self, progress, total):
        """ Updates the amount, and writes to stdout. Prints a carriage return
            first, so it will overwrite the current line in stdout."""
        print '\r',

        value = (progress * 100) / float(total)

        self.updateAmount(value)
        sys.stdout.write(str(self))
        if progress == total:
            sys.stdout.write("\n")
        sys.stdout.flush()

class IndexEntry(object):
    def __init__(self, index=None, length=None):
        self.index = index
        self.length = length
    
class FastaIndex(object):

    def __init__(self, inputFileName, rebuild=False, progress=False):

        self.indexFileName = os.path.splitext(inputFileName)[0] + ".idx"
        self.inputFile = open(inputFileName, 'r')

        # Check for existence of an index file:
        if rebuild or not os.path.exists(self.indexFileName):

            if progress:
                pg = ProgressBar()        
                fileSize = os.path.getsize(inputFileName)

            self.index = {}
            pos = 0

            self.lineWidth = None

            prevLength = 0
            prevKey = None

            for l in self._lineGenerator(self.inputFile, "\n"):

                pos += len(l)
                if progress:
                    pg.progressHook(pos, fileSize)

                if l[0] == ">":
                    key = l[1:].strip()
                    assert not self.index.has_key(key)
                    idxEnt = IndexEntry(index=pos)
                    self.index[key] = idxEnt
                    if prevKey is not None:
                        self.index[prevKey].length = prevLength

                    prevLength = 0

                    prevKey = key
                else:
                    prevLength += len(l.strip())

            if prevKey is not None:
                self.index[prevKey].length = prevLength
            
                
            indexFile = open(self.indexFileName, 'w')
            assert self.lineWidth
            pickle.dump((self.index, self.lineWidth), indexFile)
            indexFile.close()

        else:
            indexFile = open(self.indexFileName, 'r')
            self.index, self.lineWidth = pickle.load(indexFile)
            indexFile.close()

    def _lineGenerator(self, fileob, separator):
        """
        Generates lines from input file and populates the lineWidth attribute.
        """
        seplen = len(separator)
        block = fileob.read(8192)

        self.lineWidth = None

        while True:
            where = block.find(separator)

            if where < 0:
                moredata = fileob.read(8192)
                if moredata:
                    block += moredata
                    continue
                if block:
                    yield block
                return

            where += seplen

            if block[0] != ">":
                if self.lineWidth is None:
                    self.lineWidth = where-1
                assert where-1 == self.lineWidth or len(block) == where or block[where] == ">"

            yield block[:where]
            block = block[where:]


    def get(self, seq_id, start, length, strand=1):

        idxEnt = self.index[seq_id]

        assert 0 <= start < idxEnt.length and 0 < start+length <= idxEnt.length, "%d %d %d %d" % (start, idxEnt.length, start+length, idxEnt.length)

        if start == 0:
            readStart = idxEnt.index
        else:
            readStart = idxEnt.index + start + int(start/self.lineWidth)

        ## readLength = length + int(length/self.lineWidth)
        if start % self.lineWidth and length > self.lineWidth - (start % self.lineWidth):
            readLength = length + 1 + int(( length - (self.lineWidth - (start % self.lineWidth)) )/self.lineWidth)
        else:
            readLength = length + int(length/self.lineWidth)

        self.inputFile.seek(readStart)
        s = self.inputFile.read(readLength)

        s = s.replace("\n", "")

        if strand == -1:
            s = string.translate(s[::-1], string.maketrans('AGTC', 'TCAG'))

        return s

if __name__ == "__main__":

    inputFileName = sys.argv[1]

    idx = FastaIndex(inputFileName)

    if len(sys.argv) > 2:
        name = sys.argv[2]
        start = int(sys.argv[3])
        length = int(sys.argv[4])

        # snippet to print it in a nice looking way:
        width = 80
        seq = idx.get(name, start, length)          
        lines, remainder = divmod(len(seq), width)
        print ">%s_%d_%d" % (name, start, length)
        for i in range(lines):
            print seq[i*width:i*width+width]
        if remainder:
            print seq[-remainder:]
        
#     print idx.get("50262771", 0, 60)
#     print idx.get("70608782", 3, 10)
#     print idx.get("15718735", 181, 120)
