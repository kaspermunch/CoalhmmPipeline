#generic chunkenizer using 3 objects for splitting, criteria and output     
#Construct using
#Chunkenizer(
#parameter 1 : splitter - object potetially manipulating maf entries shortening or even splitting
#parameter 2 : merge criteria - a decider choosen whether or not two maf entries should be merged
#parameter 3 : output - an object responsible outputting the mafs in the chunks described by the criteria
#)
#to execute just call
#chunkenizerObject.chunkenize(
#parameter 1 : string containing the filename for the maf that should be chunkenized
#)
#WARNING: Use instances only once! depending parameter objects data may be corrupted
#TODO: consider making a use once or just a function
import sys
from MAF import *
from collections import deque
from PytablesOutput import CoordinateError

class Chunkenizer:

    def __init__(self, splitter, mergeCriteria, mafFilter, mafTest, mafQualityFilters, chunkQualityFilters, truncater, masker):
        self.mafFilter = mafFilter
        self.mafTest = mafTest
        self.splitter = splitter
        self.mergeCriteria = mergeCriteria
        self.mafQualityFilters = mafQualityFilters
        self.chunkQualityFilters = chunkQualityFilters
        self.truncater = truncater
        self.masker = masker
        self.instanceUsedAlready = False
    
    def acceptMaf(self, maf):
        for mqf in self.mafQualityFilters:
            if not mqf.accept(maf): 
                return False        
        return True
        
    def acceptChunk(self, mafs):
        for cqf in self.chunkQualityFilters:
            if not cqf.accept(mafs): 
                return False        
        return True
        
    def chunkenize(self, inputMAF, output):

        assert not self.instanceUsedAlready, "WARNING: Use instances only once! depending parameter objects data may be corrupted"
        self.instanceUsedAlready = True

        toMerge = [] #list of mafs to merge
        q = deque() #queue containing mafs which still need processing

        mafFile = open(inputMAF) 
        
        # print "DEBUGGING"
        # with open('/home/kmt/chr2.fa', 'r') as f:
        #     f.readline()
        #     chromosome2 = "".join(x.strip() for x in f.readlines())
        # from string import maketrans
        # trans = maketrans("ATGC", "TACG")

        for maf in MAFIterator(mafFile):

            self.mafFilter.inplace(maf)
                
            if not self.mafTest.test(maf):
                continue

            # print 'orig maf'
            # print maf

            #split maf and put it on the processing que
            for splitmaf in self.splitter.split(maf):
                #truncate the maf
                #print "before truncating\n", splitmaf, "\nafter"

                # print "after splitting"
                # print splitmaf
                for idx in range(maf.count()):
                    assert maf.end(idx) == maf.start(idx) + maf.length(idx)
                    assert maf.start(idx) <= splitmaf.start(idx) <= maf.start(idx) + maf.length(idx), (maf.start(idx), splitmaf.start(idx), maf.start(idx) + maf.length(idx))
                    assert maf.start(idx) <= splitmaf.end(idx) <= maf.start(idx) + maf.length(idx), (maf.start(idx), splitmaf.end(idx), maf.start(idx) + maf.length(idx))
                #     if splitmaf.name(idx).startswith('hg19'):
                #         s = splitmaf.data(idx).replace('-', '').upper()
                #         if splitmaf.strand(idx) == '-':
                #             ref = chromosome2[splitmaf.srcLength(idx) - splitmaf.start(idx) - splitmaf.length(idx): splitmaf.srcLength(idx) - splitmaf.start(idx) ].upper()
                #             ref = ref[::-1].translate(trans)
                #         else:
                #             ref = chromosome2[splitmaf.start(idx): splitmaf.start(idx) + splitmaf.length(idx)].upper()
                #         assert s == ref, (s, ref)


                splitmaf = self.truncater.truncate(splitmaf)

                #print splitmaf
                #raw_input("press enter")    

                if splitmaf != None and self.acceptMaf(splitmaf):
                    q.append(self.masker.mask(splitmaf))

                    # print "after truncation"
                    # print splitmaf
                    for idx in range(maf.count()):
                        assert maf.end(idx) == maf.start(idx) + maf.length(idx)
                        assert maf.start(idx) <= splitmaf.start(idx) <= maf.start(idx) + maf.length(idx), (maf.start(idx), splitmaf.start(idx), maf.start(idx) + maf.length(idx), idx)
                        assert maf.start(idx) <= splitmaf.end(idx) <= maf.start(idx) + maf.length(idx), (maf.start(idx), splitmaf.end(idx), maf.start(idx) + maf.length(idx), idx)
                    # print '###############'
            

            if len(q) == 0: #queue is empty get the next
                continue
            
            #if the tomerge list is empty then just add one maf from the queue
            #later code assumes toMerge is not empty
            if len(toMerge) == 0:
                toMerge.append(q.popleft())
            
            #now just flush the queue
            #for every element in the queue: test with the mergeCriteria against the last added maf in tomerge list
            #if they should merge: just add the element to the tomerge list
            #otherwise writeout the tomerge list, and then add.
            
            while len(q) > 0:
                e = q.popleft()
                last = toMerge[len(toMerge) -1]
                
                if self.mergeCriteria.shouldMerge(last, e):
                    toMerge.append(e)
                else:
                    if self.acceptChunk(toMerge):
                        try:
                            output.writeChunk(toMerge)
                        except CoordinateError, err:
                            print >>sys.stderr, err.message
                    toMerge = [e]
        
        assert output.initialized # seems all mafs were rejected        
        
        #since the tomerge list is never emty we need to write out the last chunk
        if len(toMerge) > 0 and self.acceptChunk(toMerge):
            try:
                output.writeChunk(toMerge)
            except CoordinateError, err:
                print >>sys.stderr, err.message                
        
        #fullfill the output objects last wishes... if any
        output.finalize()
