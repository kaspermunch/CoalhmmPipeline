from tables import *
from PipelineUtils import *
from IntervalQuery import *
from random import random
from math import *
from time import time as gettime

#querying pytables directly turned out to spend approx 300 ms/query
#putting the intervals in a binary tree reduced query time to micro seconds
#table parameter should be a table from hdf.root.coordinates group
#al is the alignment number
def createCoordinatesStructure(table, al):
    coordinates = IntervalQuery(notFound=(-1,-1,-1,-1))
    for x in table.where("alignmentNumber==" + str(al)):
        (_a, _b, alnumber, spbegin, spend, albegin, alend) = x.fetch_all_fields()
        coordinates.put(int(albegin), int(alend), (int(albegin), int(alend), int(spbegin), int(spend)))
        
    return coordinates

#hdf is a reference to a hdf file
#posterior table file is a string containing a path to a posterior table file
#alignmentNumber is an integer
#listNumber is an integer representing the list used to make posterior table
def importPosterior(hdf, posteriorTableFile, alignmentNumber, listNumber, coords):

    if "posteriors" not in hdf.root:
        hdf.createGroup(hdf.root, "posteriors")
    
    table = None
    if decodeChrName(alignmentNumber) in hdf.root.posteriors:
        table = hdf.root.posteriors[decodeChrName(alignmentNumber)]
    else:
        table = hdf.createTable(hdf.root.posteriors, decodeChrName(alignmentNumber), Posteriors, filters=Filters(complevel=9, complib='blosc', shuffle=True, fletcher32=False))
    data = open(posteriorTableFile)

    data.readline() #discard header
    buff = []
    
    qtime = 0.0
    starttime = gettime()
    processed = 0
    
    for x in hdf.root.lists.where("(alignmentNumber==%i)&(listNumber==%i)" % (alignmentNumber, listNumber)):
        chunk = x["chunk"]
        q = [(x2["begin"],x2["end"], x2["colId"]) for x2 in hdf.root.maps.main.where("(alignmentNumber==%i)&(chunk==%i)"% (alignmentNumber, chunk))]
        assert len(q) >= 1
        prev_end = q[0][0]
        for (alignmentStart, alignmentEnd, colid) in q:
            prev_coords = (-1,-1,-1,-1)
            for alignmentPosition in range(prev_end, alignmentEnd):
                tmp = data.readline().split(" ")
                
                states = [float(tmp[i]) for i in range(1, 5)]
                tmpstates = [(float(tmp[i]), random(),  i-1) for i in range(1, 5)]
                tmpstates.sort()
                (albegin, alend, spbegin, spend) = prev_coords
                if not(albegin <= alignmentPosition <= alend):#a small optimization
                    (albegin, alend, spbegin, spend) = coords.query(alignmentPosition)
                    prev_coords = (albegin, alend, spbegin, spend)
                    
                speciesPosition = spbegin + copysign(abs(alignmentPosition - albegin), spend-spbegin)
                
                if spend == -1:
                    speciesPosition = -1
                    
                buff.append((states[0], states[1], states[2], states[3], tmpstates[3][2], tmpstates[3][0], chunk, alignmentPosition, alignmentNumber, speciesPosition))
                
                if len(buff) >= 10000:
                    table.append(buff)
                    processed = processed + len(buff)
                    #print "processed", processed, "per sec", (processed/(gettime() - starttime))
                    del buff[:]
                    
            prev_end = alignmentEnd
            
    if len(buff) > 0:
        table.append(buff)
        processed = processed + len(buff)
#        print "processed", processed, "per sec", (processed/(gettime() - starttime))
        buff = []

""" example use
hdf = openFile("/tmp/julien.h5", "a")        
basedir = "/home/bonobo1/tmpdir/jdutheil/Data/CoalHMM/Gorilla3/CoalHMM/ILS09_gamma/Lists1Mb_aln/"

for al in range(1,23):
    nlists = len([x for x in hdf.root.lists.where("(alignmentNumber==%i)&(listIndex==1)" % al)])
    coords = createCoordinatesStructure(hdf.root.coordinates.human, al)
    for i in range(1, nlists +1):
        print "doing chr", al, "list", i, "of", nlists
        importPosterior(hdf, basedir + "/Chr" + str(al)+ "/chr" + str(al)+ "." + str(i) + ".ILS09_gamma.hmm.values.csv", al, i, coords)
"""


