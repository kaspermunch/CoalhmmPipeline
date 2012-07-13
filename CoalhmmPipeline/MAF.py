
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
                    self.alignment.append(l.split())
                elif l.startswith("q"):
                    self.quality.append(l.split())
        # sort by name
        self.alignment.sort(key=lambda x: x[1].split(".")[0])

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
    	return len(self.data(0))

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

