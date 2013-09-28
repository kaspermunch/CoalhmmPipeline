from tables import *
from PipelineUtils import *
import os.path
import os


class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass


class ListGenerationError(Error):
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


def generateLists(alignmentNumber, size, hdf):
    list_entries = []
    chunking = [x.fetch_all_fields() for x in hdf.root.maps.main.where("(alignmentNumber=="+str(alignmentNumber) +") & (segment==1)")] #assume sorted
    assert chunking, "no chunks for " + str(alignmentNumber)

    list_start = -2*size
    curr_list = 0
    list_entry = 1
    # FIXME: should we start a new list when we change chromosome?
    for (colId, alno, chunk, segment, score, begin, end) in chunking:
        if(begin - list_start > size):
            curr_list = curr_list+1
            list_start = begin
            list_entry = 1
            
            # get start and end of entire chunk not just the first segment:
            segments = [x.fetch_all_fields() for x in hdf.root.maps.main.where("(alignmentNumber==%d) & (chunk==%d)" % (alignmentNumber, chunk))]
            _, _, _, _, _, beginLst, endLst = zip(*segments)
            begin, end = min(beginLst), max(endLst)

        list_entries.append((curr_list, list_entry, alignmentNumber, chunk, begin, end))
        list_entry = list_entry +1
        
    if "lists" not in hdf.root:
        hdf.createTable(hdf.root, "lists", Lists, filters=Filters(complevel=0, complib='blosc', shuffle=True, fletcher32=False))
    elif len([x for x in hdf.root.lists.where("alignmentNumber=="+str(alignmentNumber))]) > 0:
        raise ListGenerationError(alignmentNumber, "Lists allready generated for alignment")
        #print "Lists allready generated for alignment", alignmentNumber, "Aborting"
        return

    hdf.root.lists.append(list_entries)

def writeList(alignmentNumber, listNumber, hdf, filename, chunkdir):
    if not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))
    f = open(filename, "w")
    l = [x.fetch_all_fields() for x in hdf.root.lists.where("(alignmentNumber==" + str(alignmentNumber) +") & (listNumber==" + str(listNumber) +")")]
    if not chunkdir.endswith("/"):
        chunkdir = chunkdir + "/"
    for (listNumber, listIndex, alignmentNumber, chunk, start, end) in l:
        subDir = "%d/%d/" % (alignmentNumber, int(chunk)/500)  # needs to fit same formula in PytablesOutput.py
        f.write(chunkdir + subDir + str(alignmentNumber) + "." + str(chunk) + ".fasta\n")
    f.close()
    
def numberOfLists(alignmentNumber, hdf):
    return len([x for x in hdf.root.lists.where("(alignmentNumber==" + str(alignmentNumber) + ")&(listIndex==1)")])
    
def writeAllLists(alignmentNumber, hdf, chunkdir, listdir, outputFileTemplate=''):
    if not listdir.endswith("/"):
        listdir = listdir + "/"
    #print listdir, alignmentNumber, numberOfLists(alignmentNumber, hdf)
    for l in xrange(1, numberOfLists(alignmentNumber, hdf) +1):

        if outputFileTemplate:
            # FIXME: CHANGE HERE FOR NEW RUN
            fileName = outputFileTemplate % (decodeChrName(alignmentNumber), str(l))
            #fileName = outputFileTemplate % (alignmentNumber, str(l))
        else:
            # FIXME: CHANGE HERE FOR NEW RUN
            fileName = listdir + str(decodeChrName(alignmentNumber)) + "." + str(l) + ".txt"
            #fileName = listdir + str(alignmentNumber) + "." + str(l) + ".txt"

        writeList(alignmentNumber, l, hdf, fileName, chunkdir)

# def writeAllLists(alignmentNumber, hdf, chunkdir, listdir, outputFileTemplate=''):
#     if not listdir.endswith("/"):
#         listdir = listdir + "/"
#     for l in xrange(1, numberOfLists(alignmentNumber, hdf) +1):
#         writeList(alignmentNumber, l, hdf, listdir + str(alignmentNumber) + "." + str(l) + ".txt", chunkdir)
