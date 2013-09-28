
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



class Forward:

    def __init__(self, hdf, analysis, nrStates):
        self.hdf = hdf
        self.nrStates = nrStates
        if not ("forward" in self.hdf.root):
            self.hdf.createGroup(self.hdf.root, "forward")
            self.analysis = listNr

    def _createTable(self):

        scheme = dict()
        scheme['listNumber'] = UInt16Col(pos=1)
        scheme['position'] = Int64Col(pos=2)
        for i in range(self.nrStates):
            scheme["V" + str(i)] = Float64Col(pos=(i+3))
        table = hdf.createTable(hdf.root.forward, listNr, scheme, filters=Filters(complevel=9, complib='blosc', shuffle=True, fletcher32=False))


    def _createCSIndex(self):
        if "position" not in self.hdf.root.forward[listNr].colindexes or not hdf.root.forward[listNr].colindexes["position"].is_CSI:
            hdf.root.forward[listNr].cols.position.createCSIndex()        




