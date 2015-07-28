
class MAF(object):

    def __init__(self, block=None):
        self.alignment = []
        self.quality = []
        self.annotation = ""
        if block is not None:
            for l in block.strip().split("\n"):
                if l.startswith("#"):
                    continue
                if l.startswith("a"):
                    self.annotation = l
                elif l.startswith("s"):
                    lst = l.split()
                    if len(lst) == 6:
                        # there was no sequence
                        lst.append('')
                    self.alignment.append(lst)
                elif l.startswith("q"):
                    self.quality.append(l.split())
        # sort by name
        self.alignment.sort(key=lambda x: x[1].split(".")[0])

        # check core features
        for i in range(self.count()):
            assert self.chromosome(i) is not None
            assert self.name(i) is not None
            assert self.start(i) is not None
            assert self.size(i) is not None
            assert self.strand(i) is not None
            assert self.length(i) == len(self.data(i).replace('-', '')), (self.length(i), len(self.data(i).replace('-', '')), self.name(i), self.data(i))
            assert self.size(i) == len(self.data(i)), (self.size(i), len(self.data(i)), self.name(i), self.data(i))

    def applyFilter(minScore, ignore=[]):
        """
        Applies filter on MAF inplace
        """
        assert self.quality
        for sp in self.count():
            assert len(self.quality[sp][6]) == len(self.alignment[sp][6])
            for i in self.size():
                q = self.quality[sp][6][i]
                if q not in ignore and q < minScore:
                    self.alignment[sp][6][i] = 'N'        

    def count(self):
    	return len(self.alignment)
    	
    def name(self, idx):
    	return self.alignment[idx][1].split(".")[0]

    def chromosome(self, idx):
    	tmp = self.alignment[idx][1].split(".")
    	if len(tmp) > 1:
    		c = tmp[len(tmp) -1]
    		if c.startswith("chr"):
    			return c[3:]
    		else:
    			return c
    	else:
    		return None
    
    def start(self, idx):
    	return int(self.alignment[idx][2])
    	
    def size(self, idx):
    	return len(self.data(idx))

    def length(self, idx):
    	return int(self.alignment[idx][3])
    	
    def end(self, idx):
    	return self.start(idx) + self.length(idx)
    
    def strand(self, idx):
    	return self.alignment[idx][4]
    	
    def srcLength(self, idx):
    	return int(self.alignment[idx][5])
    	
    def data(self, idx):
    	return self.alignment[idx][6]
    
    def score(self):
    	for e in self.annotation.split():
    		if e.startswith("score"):
    			return float(e.split("=")[1])
    	return -1.0
    
    def ltrunc(self, trunc_start, trunc_end):

        newmaf = MAF()

        # no meaningfull annotation
        self.annotation = "a "
        newmaf.annotation = "a "

        for i in range(len(self.alignment)):

            tag, name, start, length, strand, total, seq = self.alignment[i]

            # cutout sequences fragments:
            trunc = seq[:int(trunc_start)]
            cutout = seq[int(trunc_start):int(trunc_end)]
            self.alignment[i][6] = seq[int(trunc_end):]

            # get the number of bases we trunc off:
            offset = 0
            gaps = 0
            for c in trunc:
                if c != '-':
                    offset += 1
                else:
                    gaps += 1
                    
            # add sequence lines to new maf:
            newmaf.alignment.append([tag, name, start, str(offset), strand, total, trunc])

            # get the number of bases in the cutout:
            for c in cutout:
                if c != '-':
                    offset += 1

            # adjust coordinates on original maf:
            self.alignment[i][2] = str(int(self.alignment[i][2]) + offset)
            self.alignment[i][3] = str(int(self.alignment[i][3]) - offset)

        return newmaf

    def rtrunc(self, trunc_start, trunc_end):

        newmaf = MAF()

        # no meaningfull annotation
        self.annotation = "a "
        newmaf.annotation = "a "

        for i in range(len(self.alignment)):

            tag, name, start, length, strand, total, seq = self.alignment[i]

            # cutout sequences fragments:
            self.alignment[i][6] = seq[:int(trunc_start)]
            cutout = seq[int(trunc_start):int(trunc_end)]
            trunc = seq[int(trunc_end):]

            # get the number of bases we trunc off:
            offset = 0
            gaps = 0
            for c in trunc:
                if c != '-':
                    offset += 1
                else:
                    gaps += 1
                    
            # add sequence lines to new maf:
            newmaf.alignment.append([tag, name, start, str(offset), strand, total, trunc])

            # get the number of bases in the cutout:
            for c in cutout:
                if c != '-':
                    offset += 1

            # adjust coordinates on original maf:
            #self.alignment[i][2] = str(int(self.alignment[i][2]) + offset)
            self.alignment[i][3] = str(int(self.alignment[i][3]) - offset)

        return newmaf


    def ltrunc_by_species_coord(self, idx, sp_start, sp_end):
        """Trim maf using plus strand coordinates of a species"""

        newmaf = MAF()

        # no meaningfull annotation
        self.annotation = "a "
        newmaf.annotation = "a "

        tag, name, start, length, strand, total, seq = self.alignment[idx]        

        if strand == '-':
            start, end = int(total) - (int(start)+int(length)), int(total) - int(start)
            seq = seq[::-1]
        else:
            start, end = int(start), int(start)+int(length)

        assert start < end
        assert start <= sp_start < end
        assert start < sp_end <= end

        # number of bases to trim off start:
        trim_start = sp_start - start
        trim_end = end - sp_end

        cols_trim_start = 0
        bases = 0
        for c in seq:
            if bases == trim_start:
                break
            if c != '-':
                bases += 1
            cols_trim_start += 1
        assert bases ==  trim_start

        # number of bases to trim off end:
        cols_trim_end = 0
        bases = 0
        for c in seq:
            if bases == trim_end:
                break
            if c != '-':
                bases += 1
            cols_trim_end += 1
        assert bases ==  trim_end

        if strand == '-':
            cols_trim_start, cols_trim_end = cols_trim_end, cols_trim_start
            return self.rtrunc(cols_trim_start, len(seq)-cols_trim_end)
        else:
            return self.ltrunc(cols_trim_start, len(seq)-cols_trim_end)


    def prepend(self, other):
        for i in range(len(self.alignment)):
            assert self.alignment[i][0] == other.alignment[i][0] and self.alignment[i][1] == other.alignment[i][1]
            self.alignment[i][6] = other.alignment[i][6] + self.alignment[i][6]
            self.alignment[i][2] = self.alignment[i][2] - other.alignment[i][3]
            self.alignment[i][3] = self.alignment[i][3] + other.alignment[i][3]
        
    def __str__(self):
        s = self.annotation + "\n"
        s += "\n".join([" ".join(x) for x in self.alignment]) + "\n"
        return s + "\n\n"
    

def MAFIterator(fileob):

    openedFile = False
    if not isinstance(fileob, file):
        fileob = open(fileob)
        openedFile = True

    separator = "\n\n" #unix lineendings only
    seplen = len(separator)
    block = fileob.read(8192)

    while True:
        where = block.find(separator)

        if where < 0:
            moredata = fileob.read(8192)
            if moredata:
                block += moredata
                continue
            if block:
                yield MAF(block)
            if openedFile:
                fileob.close()
            return

        where += seplen
        yield MAF(block[:where])
        block = block[where:]

if __name__ == "__main__":

    test_maf = """a score=15970.000000
s calJac3.chr7 39980538 10 + 155834243 GCGGTGTTGT
s gorGor3.chr1 5749428 10 + 229507203 GTGGCGTTGC
s hg19.chr1 5682968 10 + 249250621 GTGGCGTTGC
s nomLeu1.GL397433 2526565 10 + 2715696 GTGGCGTTGT
s panTro3.chr1 5590671 10 + 228333871 GTGGCGTTGC
s ponAbe2.chr1 5012978 10 - 229942017 GTGGCATTGT
s rheMac2.chr1 8619867 10 + 228252215 GTGGTGTTGT"""

    maf = MAF(test_maf)
    print maf
    #print maf.ltrunc_by_species_coord(2, 5682968+2, 5682978-3)
    print maf.ltrunc_by_species_coord(2, 5682968, 5682968+4)
    print maf

    print "########"

    test_maf = """a score=4418638.000000
s calJac3.chr7 115851821 1872 - 155834243 ttttctagtttcttaaggtggaagctgaagttgttgatGGCACtcagtccaagactgattattccccatgattgaagc
s gorGor3.chr1 5722038 1873 + 229507203 ttttctagtttcttaaggtggaagctgaagttgttgatGGCACtcagtccaagactgattattccccatgattgaagc
s hg19.chr1 243565767 68 - 249250621 ttttctagtttcttaaggtggaagctgaagttgttgatGGCACtcagtccaggactga----------tgattgaagc
s nomLeu1.GL397433 187229 1891 - 2715696 ttttctagtttcttaaggtggaagctgaagttgttgaCGGCACtcagtccaggactgattattccccatgattgaggc
s panTro3.chr1 222741314 1875 - 228333871 ttttctagtttcttaaggtggaagctgaagttgttgaCGGCACtcagtccaggactgattattccccatgactgaagc
s papHam1.scaffold1958 6904 1960 - 191851 ttttctagtttcttaaggtggaagctgaagtcgttgaCGTCACtcagtccaggactaattattccccatgactgaagc
s ponAbe2.chr1 224927154 1874 + 229942017 ttttctagtttcttaaggtggaagctgaagttgttgatGGCACtcagtccaggactgattattccccatgactgaagc
s rheMac2.chr1 219630387 1949 - 228252215 ttttctagtttcttaaggtggaagctgaagttgttgaCGCCACtcagtccaggactaattattccccatgactgaagc"""

    maf = MAF(test_maf)
    print maf
    
    #print maf.ltrunc_by_species_coord(2, 5684787+20, 5684855-10)
    print maf.ltrunc_by_species_coord(2, 5684787, 5684787+4)
    print maf


