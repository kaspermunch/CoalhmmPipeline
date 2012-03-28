from tables import *
from PipelineUtils import *

def fillGaps(table, alignmentNumber):
    prev_end = -1
    buff = []
    for x in table.itersorted(table.cols.alignmentPositionBegin):
        (chromosome, strand, alno, spbegin, spend, albegin, alend) = x.fetch_all_fields()
        if alno != alignmentNumber:
            continue
            
        if prev_end +1 == albegin:
            prev_end = alend
            continue
            
        buff.append((chromosome, strand, alno, -1,-1, prev_end +1, albegin-1))
        prev_end = alend
        
        if len(buff) % 1000 == 0:
            print "inserted", len(buff)
        
    buff.append((0, 0, 0, -1,-1, prev_end +1, prev_end+0xFFFFFFFFFFFF))
    return buff
    
    
hdf = openFile("/tmp/test3.h5", "a")
buffs = []
for i  in range(1,24):
    buffs.append(fillGaps(hdf.root.coordinates.human, i))
    print "done", i
    
for x in buffs:
    hdf.root.coordinates.human.append(x)

print "reindexing albegin"
hdf.root.coordinates.human.cols.alignmentPositionBegin.reIndex()
print "reindexing alend"
hdf.root.coordinates.human.cols.alignmentPositionEnd.reIndex()
print "reindexing alno"
hdf.root.coordinates.human.cols.alignmentNumber.reIndex()

hdf.close()
