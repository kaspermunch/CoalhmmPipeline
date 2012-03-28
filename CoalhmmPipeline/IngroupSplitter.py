from MAF import MAF

class IngroupSplitter:
    def __init__(self, ingroup, splitdist):
        self.ingroup = ingroup
        self.splitdist = splitdist
    
    #return a list of mafs representing the old mafs that have split
    #or a list containing just the original maf if no splitting is required    
    def split(self, maf):
        in_al = []
        in_pos = []
        index = []
        indir = []
        out_al = []
        out_pos = []
        outdex = []
        outdir = []
        for i in range(maf.count()):
            if maf.name(i) in self.ingroup:
                in_al.append(maf.data(i))
                in_pos.append(maf.start(i))
                index.append(i)
                if maf.strand(i) == '+':
                    indir.append(1)
                else:
                    indir.append(-1)
            else:
                out_al.append(maf.data(i))
                out_pos.append(maf.start(i))
                outdex.append(i)
                if maf.strand(i) == '+':
                    outdir.append(1)
                else:
                    outdir.append(-1)
                
        gapcount = 0
        for i in range(len(in_al[0])):
            ns = 0
            
            for j in range(len(in_al)):
                if in_al[j][i] not in "-":
                    in_pos[j] = in_pos[j] +1
                
                if in_al[j][i] in 'N-':
                    ns = ns + indir[j]
                    
            for j in range(len(out_al)):
                if out_al[j][i] not in "-":
                    out_pos[j] = out_pos[j] + outdir[j]
                    
            if ns == len(in_al):
                gapcount = gapcount +1
            elif gapcount <= self.splitdist:
                gapcount = 0
            elif gapcount > self.splitdist:
                maf1 = "a score=%f\n" % maf.score()
                for x in range(len(index)):
                    ind = index[x]
                    start = maf.start(ind)
                    size = abs(in_pos[x] - start) - gapcount
                    if in_al[x][i] not in "-":
                        size = size -1
                    maf1 += "s %s.%s %i %i %s %i %s\n" % (maf.name(ind), maf.chromosome(ind), start, size, maf.strand(ind), maf.srcLength(ind), maf.data(ind)[0:i-gapcount])
                    
                for x in range(len(outdex)):
                    ind = outdex[x]
                    start = maf.start(ind)
                    size = abs(out_pos[x] - start) - gapcount
                    if out_al[x][i] not in "-":
                        size = size -1
                    maf1 += "s %s.%s %i %i %s %i %s\n" % (maf.name(ind), maf.chromosome(ind), start, size, maf.strand(ind), maf.srcLength(ind), maf.data(ind)[0:i-gapcount])
                
                maf2 = "a score=%f\n" % maf.score()
                for x in range(len(index)):
                    ind = index[x]
                    start = in_pos[x]
                    if in_al[x][i] not in "-":
                        start = start - indir[x]
                        
                    size = maf.length(ind) - (start - maf.start(ind))
                    maf2 += "s %s.%s %i %i %s %i %s\n" % (maf.name(ind), maf.chromosome(ind), start, size, maf.strand(ind), maf.srcLength(ind), maf.data(ind)[i:])
                    
                for x in range(len(outdex)):
                    ind = outdex[x]
                    start = out_pos[x]
                    if out_al[x][i] not in "-":
                        start = start - outdir[x]
                    size = maf.length(ind) - (start - maf.start(ind))
                    maf2 += "s %s.%s %i %i %s %i %s\n" % (maf.name(ind), maf.chromosome(ind), start, size, maf.strand(ind), maf.srcLength(ind), maf.data(ind)[i:])    
                gapcount = 0
                
                tmpmaf = MAF(maf1)
                res = []
                
                if(len(tmpmaf.data(0)) > 0):
                    res.append(tmpmaf)
                
                for m in self.split(MAF(maf2)):
                    res.append(m)
                    
                return res
        
        return [maf]
