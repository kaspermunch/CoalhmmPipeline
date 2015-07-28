from MAF import MAF

class IngroupSplitter:
    def __init__(self, ingroup, splitdist, junkchars='Nn'):
        self.ingroup = ingroup
        self.splitdist = splitdist
        self.junkchars = junkchars
        assert "-" not in self.junkchars

    
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

            assert maf.start(i) <= maf.end(i), (maf.start(i), maf.end(i))

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

        # keep track of how many gaps there are in each sequence in current block of junk
        block_dashes = dict()
        for j in index + outdex:
            block_dashes[j] = 0

        for i in range(len(in_al[0])):
            ns = 0

            for j in range(len(in_al)):
                if in_al[j][i] not in "-":
                    in_pos[j] = in_pos[j] +1
                else:
                    block_dashes[index[j]] += 1
                
                if in_al[j][i] in self.junkchars + '-':
                    ns = ns + 1#indir[j]

            for j in range(len(out_al)):
                if out_al[j][i] not in "-":
                    out_pos[j] = out_pos[j] + 1#outdir[j]
                else:
                    block_dashes[outdex[j]] += 1
                    
            if ns == len(in_al):
                gapcount = gapcount +1
            elif gapcount <= self.splitdist:
                gapcount = 0
                for j in index + outdex:
                    block_dashes[j] = 0
            elif gapcount > self.splitdist:
                maf1 = "a score=%f\n" % maf.score()
                for x in range(len(index)):
                    ind = index[x]
                    start = maf.start(ind)
                    # size = abs(in_pos[x] - start) - gapcount
                    assert in_pos[x] >= start
#                    assert gapcount >= block_dashes[ind], (gapcount, block_dashes[ind], ind)
                    size = in_pos[x] - start - (gapcount - block_dashes[ind])                    
                    # if in_al[x][i] not in "-":
                    size = size -1
                    # print "#### in:", in_pos[x], gapcount, block_dashes[ind], start, size                        
                    maf1 += "s %s.%s %i %i %s %i %s\n" % (maf.name(ind), maf.chromosome(ind), start, size, maf.strand(ind), maf.srcLength(ind), maf.data(ind)[0:i-gapcount])
                    
                for x in range(len(outdex)):
                    ind = outdex[x]
                    start = maf.start(ind)
                    size = out_pos[x] - start - (gapcount - block_dashes[ind])
                    # if out_al[x][i] not in "-":
                    size = size -1
                    # print "#### out:", out_pos[x], gapcount, block_dashes[ind], start, size
                    maf1 += "s %s.%s %i %i %s %i %s\n" % (maf.name(ind), maf.chromosome(ind), start, size, maf.strand(ind), maf.srcLength(ind), maf.data(ind)[0:i-gapcount])
                
                maf2 = "a score=%f\n" % maf.score()
                for x in range(len(index)):
                    ind = index[x]
                    start = in_pos[x]
                    if in_al[x][i] not in "-":
                        start = start - 1#indir[x]
                    size = maf.length(ind) - (start - maf.start(ind))
                    maf2 += "s %s.%s %i %i %s %i %s\n" % (maf.name(ind), maf.chromosome(ind), start, size, maf.strand(ind), maf.srcLength(ind), maf.data(ind)[i:])
                    
                for x in range(len(outdex)):
                    ind = outdex[x]
                    start = out_pos[x]
                    if out_al[x][i] not in "-":
                        start = start - 1#outdir[x]
                    size = maf.length(ind) - (start - maf.start(ind))
                    maf2 += "s %s.%s %i %i %s %i %s\n" % (maf.name(ind), maf.chromosome(ind), start, size, maf.strand(ind), maf.srcLength(ind), maf.data(ind)[i:])    
                gapcount = 0

                tmpmaf = MAF(maf1)
                res = []
                # print "## 1 ##"
                # print tmpmaf 
                for x in range(tmpmaf.count()):
                    assert tmpmaf.start(x) <= tmpmaf.start(x), (tmpmaf.start(x), tmpmaf.end(x))

                if(len(tmpmaf.data(0)) > 0):
                    res.append(tmpmaf)

                tmpmaf2 = MAF(maf2)
                # print "## 2 ##"
                # print tmpmaf2                
                for x in range(tmpmaf2.count()):
                    assert tmpmaf2.start(x) <= tmpmaf2.start(x), (tmpmaf2.start(x), tmpmaf2.end(x))

                if(len(tmpmaf2.data(0)) > 0):
                    for m in self.split(tmpmaf2):
                        res.append(m)
#                 for m in self.split(MAF(maf2)):
#                     res.append(m)
                    
                return res
        
        return [maf]
