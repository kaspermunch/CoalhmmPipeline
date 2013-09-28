#Output to pytables file
#object constructed using
#PytablesOutput(
#parameter 1 : string with a path to a directory where the fasta files should be stored
#parameter 2 : list containing names of ingroups ie. ["hg18", bonobo"]
#parameter 3 : string with a filename for the pytables file
#)

from tables import *
import os
import sys
import glob
from MAF import *
from PipelineUtils import *

class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass


class CoordinateError(Error):
    """
    Exception raised for errors in ...

    Attributes:
        expression: input expression in which
                    the error occurred
        message:    explanation of the error
    """
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


#class responsible for the entire output
#this is where nearly all the work goes
class PytablesOutput:
    def __init__(self, chunkDir, ingroup, tableHDF, alignmentNumber, chunkFilePrefix=''):
        self.initialized = False
        self.nextChunk = 1 #internal counter to keep track if chunk numbering
        self.chunkDir = chunkDir #output directory for the chunked fasta files
        self.segment = 1 #counts which maf we are dealing with in a chunk. Is frequntly reset
        self.currentPosition = 0 #the current position in the alignment
        self.ingroup = ingroup #list of the ingroup species
        self.ingroupIndexes = None #speedup technique. All species in the mafs have indexes form 0..n-1, here we keep the indexes of those
                                   #that are in the ingroup
        self.coordinateTables = None #stores references to tables containing coordinates
                                     #they are stored at the same index as in the maf
        self.alignmentNo = alignmentNumber #alignment number stored together with the data, typically follows chromosome
        self.mafCount = 0 #running counter. Counts how many mafs have currently been processed.
                          #it is also used to give unique row ids to the maps

        self.hdfFile = tableHDF
        #Create two groups
        if not ("maps" in self.hdfFile.root):
            self.hdfFile.createGroup(self.hdfFile.root, "maps")
        if not ("coordinates" in self.hdfFile.root):
            self.hdfFile.createGroup(self.hdfFile.root, "coordinates")
        
        #we know we are going to need this table, so we will create it now
        if not ("main" in self.hdfFile.root.maps):
            self.hdfFile.createTable(self.hdfFile.root.maps, "main", MainChunkMap)

    #We don't know what the mafs will look like until we get the first one
    #here we initialize some structures to help output the result
    def initialize(self, maf):
        #speed up trick, replaces something like "name in ingroup"
        #with just traversing this list of integers
        #speeds up alot
        self.ingroupIndexes = []
        for i in range(maf.count()):
            if maf.name(i) in self.ingroup:
                self.ingroupIndexes.append(i)

        
        #it takes a long time to look up tables by their name
        #so we store a reference in a list    
        self.coordinateTables = []
        for i in range(maf.count()):
            self.hdfFile.createTable(self.hdfFile.root.maps, maf.name(i), SpeciesChunkMap)
            if i not in self.ingroupIndexes: #outgroup species do not get a table
                self.coordinateTables.append(None) #but outgroups still have an index we should respect
            else:
                #create a table for coordinates
                self.coordinateTables.append([]) # rows gets inserted as tubles and flushed to pytables later
                if not (maf.name(i) in self.hdfFile.root.coordinates):
                    self.hdfFile.createTable(self.hdfFile.root.coordinates, maf.name(i), SpeciesCoordinates, filters=Filters(complevel=9, complib='blosc', shuffle=True, fletcher32=False)) #creation of the table

    
    #this method gets called with the mafs (in a list) that should be chunked together
    def writeChunk(self, mafs):
        #test if this is the first mafs we receive for output
        #if so we we should initialize the rest of the structures
        if not self.initialized:
            self.initialize(mafs[0])
            self.initialized = True

        #report progress
#         print "merging", len(mafs), "mafs, total processed", self.mafCount
    
        self.segment = 1 #reset segment counter! we are starting on a new chunk!

#         for m in mafs:
#             print 'PytablesOutput, chunk len:', m.size(0)

        self.writeFasta(mafs) #todo: consider inlining this... it's 30 lines of code

        self.nextChunk = self.nextChunk +1 
        
    def finalize(self):
        #print "Creating indexes"
        #for name in self.ingroup:
        #    print name
            #self.hdfFile.root.coordinates._f_getChild(name).cols.alignmentPosition.createIndex()
            #self.hdfFile.root.coordinates._f_getChild(name).cols.SpeciesPositionOnPlusStrand.createIndex()
        #self.hdfFile.close()
        return
            
    def writeFasta(self, mafs):
        prevMAF = None #we are merging the first one
        data = [] #Each entry is this list corresponds to a species. A species is in turn represented by a list of data
        dataKeys = [] #name of the species in the data list. indexes matched
        #we the initialize these variables
        for i in range(mafs[0].count()):
            data.append([])
            dataKeys.append(mafs[0].name(i))
            
        for maf in mafs:
            if prevMAF == None: #is this the first maf
                #then just append the data and move on to the next maf
                prevMAF = maf
                self.appendToChunk(data, maf)
                continue

            ##insert gap between prevMAF and maf in output
            for i in range(maf.count()):
                data[i].extend(self.getGap(prevMAF, maf, i)) # FIXME: getGap takes max gap not gap for i (fine if only if they are always the same...)

# #######################################################
#                 # FIXME: here we should add the gap to the beginning of the maf before appending to chunk
#                 gap = self.getGap(prevMAF, maf, i)
#                 mafstr = ''
#                 for j in  range(maf.count()):
#                     mafstr += "s %s.%s %i %i %s %i %s\n" % (maf.name(j), maf.chromosome(j), maf.end(j), len(gap), maf.strand(j), maf.srcLength(j), gap)
#                 maf.prepend(MAF(mafstr))
# #######################################################
                
            
            self.currentPosition += len(self.getGap(prevMAF, maf, 0)) #update the alignment position according to this gap
            
            #append maf
            self.appendToChunk(data, maf)
                    
            prevMAF = maf

        #write chunk as fasta
        subdirName = "%d/%d" % (self.alignmentNo, self.nextChunk/500) # needs to fit same formula in GenerateLists.py
        chunkDir = os.path.join(self.chunkDir, subdirName)
        chunkFile = self.openFile(chunkDir, str(self.alignmentNo) + "."+ str(self.nextChunk) + ".fasta") #open a file to write out the fasta
#         #write chunk as fasta
#         chunkFile = self.openFile(self.chunkDir, str(self.alignmentNo) + "."+ str(self.nextChunk) + ".fasta") #open a file to write out the fasta

        for i in range(len(data)):
            chunkFile.write(">"+dataKeys[i]+"\n")
            chunkFile.write("".join(data[i]) + "\n\n") #todo: consider if it would be faster to skip the join
        chunkFile.close()
    
    def appendToChunk(self, data, maf):
        #add to table where in the alignment we start
        #these data are used in main map table
        alignmentMap = (self.mafCount, self.alignmentNo, self.nextChunk, self.segment, maf.score(), self.currentPosition)
        
        self.segment += 1 #increment segment counter. we increment here because of the use in the above line
                          #it in not used anywhere else except when resetting it
        
        
        speciesPosition = [] #used in the coordinate tables to keep track of the current position in indexed species
        for i in range(maf.count()):
                speciesPosition.append(maf.start(i))
    
        
        #run over every position in the alignment to write the coordinate tables
        #and append data to the fasta
        for pos in xrange(len(maf.data(0))):
            #test if all ingroup species have a gap at position 'pos'
            #if they have, skip this position
            gaps = 0
            for i in self.ingroupIndexes:
                if maf.data(i)[pos] == "-":
                    gaps += 1
                    
            if gaps == len(self.ingroup): #skip this position, all ingroup species have a gap
                continue
                
            #append this position to the coordinate tables
            for i in self.ingroupIndexes:

                if maf.data(i)[pos] == "-": #gap/nothing does not have a position in a species, so skip
                    continue
                    
                if maf.strand(i) == "-": #all coordinates in the table should be on the plus strand
                    speciesPositionPlusStrand = maf.srcLength(i) - speciesPosition[i]
                else:
                    speciesPositionPlusStrand = speciesPosition[i]

                if speciesPositionPlusStrand < 0:
                   raise CoordinateError(speciesPositionPlusStrand < 0, "logically inconsistent input:\n%s\n" % '\n'.join([x[:80] for x in str(maf).split('\n')]))
                    
#                     # had to condition the exception because only Hsap coordinates are checked in Kays bonobo alignments. 
#                     if maf.name(i) == 'Hsap':
#                         raise CoordinateError(speciesPositionPlusStrand < 0, "logically inconsistent input: sourcelength: %d position: %d"  % (maf.srcLength(i), speciesPosition[i]))


#                     # FIXME: had to take the exception out to allow juliens bug to pass for the non-refs
#                     if maf.name(i) == 'Hsap':
#                         print "logically inconsistent input for %s: sourcelength: %d position: %d"  % (name.name(i), maf.srcLength(i), speciesPosition[i])
#                     else:
#                         raise CoordinateError(speciesPositionPlusStrand < 0, "logically inconsistent input: sourcelength: %d position: %d"  % (maf.srcLength(i), speciesPosition[i]))
                
                speciesPositionPlusStrand = speciesPositionPlusStrand & 0xFFFFFFFFFFFF #clamp to make fit in a 64 bit variable
                
                #int casts only in case of errors in the input    
                intstrand = 0
                if maf.strand(i) == "+":
                    intstrand = 1
                else:
                    intstrand = -1
                self.coordinateTables[i].append((encodeChrName(maf.chromosome(i)),  intstrand, int(speciesPositionPlusStrand), self.alignmentNo, int(self.currentPosition)))    #append to coordinate table
                speciesPosition[i] += 1
            
            #append data to the fasta file                    
            for i in range(maf.count()):    
                data[i].append( maf.data(i)[pos] )
                
            #update our current position in the alignment
            self.currentPosition += 1
        
        #flush coordinate table to the database file
        for i in self.ingroupIndexes:
            self.hdfFile.root.coordinates[maf.name(i)].append(self.coordinateTables[i])
            del self.coordinateTables[i][:]
        
        #we are missing the end component in the alignment/main map because we didn't know it
        #we can now add it, and write it to the table
        (_a,_b,_c,_d,_e, _f) = alignmentMap #_a,_b... just used for shortening the code
        alignmentMap = (_a,_b,_c,_d,_e, _f, self.currentPosition)
        self.hdfFile.root.maps.main.append([alignmentMap])
        
        #and the rest of the maps info from the species...
        for i in range(maf.count()):
                intstrand = 0
                if maf.strand(i) == "+":
                    intstrand = 1
                else:
                    intstrand = -1
                speciesMap = (self.mafCount, encodeChrName(maf.chromosome(i)), intstrand, maf.start(i) & 0xFFFFFFFFFFFF, maf.end(i) & 0xFFFFFFFFFFFF)
                #self.hdfFile.root.maps[maf.name(i)].append(speciesMap)
                table = self.hdfFile.root.maps[maf.name(i)]
                table.append([speciesMap])
                
        self.mafCount = self.mafCount +1
    
    #calculates a gap represented as a string between to mafs
    #the interface supports that it could work for a particular species
    #but right now it just does the same
    def getGap(self, prevMAF, maf, idx):
        gap = 0
        for i in range(maf.count()):
            if maf.name(i) not in self.ingroup:
                continue
            gap = max(gap, maf.start(i) - prevMAF.end(i))
        
        return "N"*gap
    
    #opens a file given a directory and a filename    
    def openFile(self, dir, name, options="w"):
        if dir != "" and dir[len(dir) -1] != "/":
            dir = dir + "/"
            
        if dir != "" and not os.path.exists(dir):
            os.makedirs(dir)
        
        return open(dir + name, options)
