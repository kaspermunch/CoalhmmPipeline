from tables import *
##monkey patching pytables
Node.__getitem__ = lambda self, key : self._f_getChild(key)
##end of patch

def encodeChrName(name):
    if name[0:3].lower() == "chr":
	name = name[3:]
    if name.isdigit():
	return int(name)
    if len(name) == 1:
	return ord(name.upper())
    if len(name) == 2:
	return ord(name[1].upper()) * 256 + ord(name[0].upper())
    return 0

def decodeChrName(number):
    if number < ord('A'):
	return 'chr' + str(number)
    if number < 256: # one letter
        return 'chr' + chr(number)
    return 'chr' + chr(number % 256) + chr(number / 256)


#table scheme for main map table
class MainChunkMap(IsDescription):
    colId = Int64Col(pos=1)
    alignmentNumber = UInt16Col(pos=2)
    chunk = UInt32Col(pos=3)
    segment = UInt32Col(pos=4)
    score = Float64Col(pos=5)
    begin = Int64Col(pos=6)
    end = Int64Col(pos=7)

#table scheme for the species mapping
class SpeciesChunkMap(IsDescription):
    colId = Int64Col(pos=1)
    chromosome = UInt16Col(pos=2)
    strand = Int8Col(pos=2)
    begin = Int64Col(pos=4)
    end = Int64Col(pos=5)

#table scheme for the coordinates    
class SpeciesCoordinates(IsDescription):
    chromosome = UInt16Col(pos=1)
    strand = Int8Col(pos=2)
    SpeciesPositionOnPlusStrand = Int64Col(pos=3)
    alignmentNumber = UInt16Col(pos=4)
    alignmentPosition = Int64Col(pos=5)
    
#table scheme for the coordinates    
class SpeciesIntervalCoordinates(IsDescription):
    chromosome = UInt16Col(pos=1)
    strand = Int8Col(pos=2)
    alignmentNumber = UInt16Col(pos=3)
    SpeciesPositionOnPlusStrandBegin = Int64Col(pos=4)
    SpeciesPositionOnPlusStrandEnd = Int64Col(pos=5)
    alignmentPositionBegin = Int64Col(pos=6)
    alignmentPositionEnd = Int64Col(pos=7)
    
#Table scheme for storing posterior probs
class Posteriors(IsDescription):
    V0 = Float64Col(pos=1)
    V1 = Float64Col(pos=2)
    V2 = Float64Col(pos=3)
    V3 = Float64Col(pos=4)
    maxstate = UInt16Col(pos=5)
    maxP = Float64Col(pos=6)
    chunk = UInt32Col(pos=7)
    alignmentPosition = Int64Col(pos=8)
    alignmentNumber = UInt16Col(pos=9)
    speciesPostion = Int64Col(pos=10)
    
#table scheme for lists of chunks
class Lists(IsDescription):
    listNumber = UInt32Col(pos=1)
    listIndex = UInt32Col(pos=2)
    alignmentNumber = UInt16Col(pos=3)
    chunk = UInt32Col(pos=4)

#table scheme for segments
class Segments(IsDescription):
    AlignmentPositionFrom = Int64Col(pos=1)
    AlignmentPositionTo = Int64Col(pos=2)
    SpeciesPositionFrom = Int64Col(pos=3)
    SpeciesPositionTo = Int64Col(pos=4)
    State = Int8Col(pos=5)

