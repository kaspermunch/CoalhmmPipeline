from tables import *
##monkey patching pytables
Node.__getitem__ = lambda self, key : self._f_get_child(key)
##end of patch

def encodeChrName(name):
    if name[0:3].lower() == "chr":
        name = name[3:]

    if name.isdigit():
        enc = int(name)
        assert enc < 1000000 # not to interfere with contig encodings
        return enc
    elif name.startswith('C') and name[1:].isdigit():
        # contig e.g. C99999
        return int(name[1:]) + 1000000
    else:
        enc = sum(ord(c) * 256**x for c, x in zip(name, range(len(name))))
        assert enc < 2**32, (name, enc)  # needs to fit in a 32 bit signed long
        return -enc # return neg value


def decodeChrName(number):
    if number < 0:
        x = -number     
        decoded = ''
        while x:
            decoded += chr(x % 256)
            x /= 256
        return 'chr{}'.format(decoded)
    elif number > 1000000:
        # contig
        return 'C{}'.format(number - 1000000)
    else:
        # chromosome
        return 'chr{}'.format(number)

# def encodeChrName(name):
#     assert name is not None 
#     #if name[0:3].lower() in ("chr", "scf"):
#     if name[0:3].lower() == "chr":
# 	name = name[3:]
#     if name.isdigit():
# 	return int(name)
#     if len(name) == 1:
# 	#return ord(name.upper())
#         return ord(name)
#     if len(name) == 2:
# 	#return ord(name[1].upper()) * 256 + ord(name[0].upper())
#         return ord(name[1]) * 256 + ord(name[0])
#     return 0

# def decodeChrName(number):
#     if number < ord('A'):
# 	return 'chr' + str(number)
#     if number < 256: # one letter
#         return 'chr' + chr(number)
#     if number < ord('z') * 256 + ord('z'): # a chrom name
#         return 'chr' + chr(number % 256) + chr(number / 256)
#     else:
#         return 'scf' + str(number)

#table scheme for main map table
class MainChunkMap(IsDescription): # maps.main
    colId = Int64Col(pos=1) # key
    alignmentNumber = Int32Col(pos=2) # alignment nr (encoded chr nr)
    chunk = UInt32Col(pos=3) # chunk nr
    segment = UInt32Col(pos=4) # maf nr nr in chunk (for tracking chunkinization procedure)
    score = Float64Col(pos=5) # maf alignment score
    begin = Int64Col(pos=6) # alignment coord
    end = Int64Col(pos=7)   # alignment coord

#table scheme for the species mapping
class SpeciesChunkMap(IsDescription): # maps.<speciesname>
    colId = Int64Col(pos=1) # key
    chromosome = Int32Col(pos=2) # species chromosome
    strand = Int8Col(pos=3) # ect..
    begin = Int64Col(pos=4)
    end = Int64Col(pos=5)
    srcLength = Int64Col(pos=6)

#table scheme for the coordinates    
class SpeciesCoordinates(IsDescription): # coords.<speciesname>  (single pos coordinates, mapping between alignment and species coordinates)
    chromosome = Int32Col(pos=1)
    strand = Int8Col(pos=2)
    SpeciesPositionOnPlusStrand = Int64Col(pos=3)
    alignmentNumber = Int32Col(pos=4)
    alignmentPosition = Int64Col(pos=5)
    
#table scheme for the coordinates    
class SpeciesIntervalCoordinates(IsDescription): # coords.<speciesname> same as above but for intervals
    chromosome = Int32Col(pos=1)
    strand = Int8Col(pos=2)
    alignmentNumber = Int32Col(pos=3)
    SpeciesPositionOnPlusStrandBegin = Int64Col(pos=4)
    SpeciesPositionOnPlusStrandEnd = Int64Col(pos=5)
    alignmentPositionBegin = Int64Col(pos=6)
    alignmentPositionEnd = Int64Col(pos=7)
    
#Table scheme for storing posterior probs
class Posteriors(IsDescription): # posteriors.<speciesname> 
    V0 = Float64Col(pos=1) # state 0
    V1 = Float64Col(pos=2) # etc...
    V2 = Float64Col(pos=3)
    V3 = Float64Col(pos=4)
    maxstate = UInt16Col(pos=5) # state max prob
    maxP = Float64Col(pos=6)  # map prob
    chunk = UInt32Col(pos=7) # chunk this originates from
    alignmentPosition = Int64Col(pos=8) # (alignment pos)
    alignmentNumber = Int32Col(pos=9)  # alignment nr
    speciesPosition = Int64Col(pos=10)  # mapped species position 
    
#table scheme for lists of chunks
class Lists(IsDescription):
    listNumber = UInt32Col(pos=1)
    listIndex = UInt32Col(pos=2)
    alignmentNumber = Int32Col(pos=3)
    chunk = UInt32Col(pos=4)
    alignmentPositionBegin = Int64Col(pos=5)
    alignmentPositionEnd = Int64Col(pos=6)


#table scheme for segments
class Segments(IsDescription):
    AlignmentPositionFrom = Int64Col(pos=1)
    AlignmentPositionTo = Int64Col(pos=2)
    SpeciesPositionFrom = Int64Col(pos=3)
    SpeciesPositionTo = Int64Col(pos=4)
    State = Int8Col(pos=5)

