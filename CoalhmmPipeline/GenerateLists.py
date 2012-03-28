from tables import *
from PipelineUtils import *
import os.path
import os


def generateLists(alignmentNumber, size, hdf):
    list_entries = []
    chunking = [x.fetch_all_fields() for x in hdf.root.maps.main.where("(alignmentNumber=="+str(alignmentNumber) +") & (segment==1)")] #assume sorted
    list_start = -2*size
    curr_list = 0
    list_entry = 1
    for (colId, alno, chunk, segment, score, begin, end) in chunking:
        if(begin - list_start > size):
            curr_list = curr_list+1
            list_start = begin
            list_entry = 1
        list_entries.append((curr_list, list_entry, alignmentNumber, chunk))
        list_entry = list_entry +1
        
    if "lists" not in hdf.root:
#        print "creating lists table"
        hdf.createTable(hdf.root, "lists", Lists, filters=Filters(complevel=0, complib='blosc', shuffle=True, fletcher32=False))
    elif len([x for x in hdf.root.lists.where("alignmentNumber=="+str(alignmentNumber))]) > 0:
        print "Lists allready generated for alignment", alignmentNumber, "Aborting"
        return
            
    hdf.root.lists.append(list_entries)

def writeList(alignmentNumber, listNumber, hdf, filename, chunkdir):
    if not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))
    f = open(filename, "w")
    l = [x.fetch_all_fields() for x in hdf.root.lists.where("(alignmentNumber==" + str(alignmentNumber) +") & (listNumber==" + str(listNumber) +")")]
    if not chunkdir.endswith("/"):
        chunkdir = chunkdir + "/"
    for (listNumber, listIndex, alignmentNumber, chunk) in l:
        f.write(chunkdir + str(alignmentNumber) + "." + str(chunk) + ".fasta\n")
    f.close()
    
def numberOfLists(alignmentNumber, hdf):
    return len([x for x in hdf.root.lists.where("(alignmentNumber==" + str(alignmentNumber) + ")&(listIndex==1)")])
    
def writeAllLists(alignmentNumber, hdf, chunkdir, listdir, outputFileTemplate=''):
    if not listdir.endswith("/"):
        listdir = listdir + "/"
    for l in xrange(1, numberOfLists(alignmentNumber, hdf) +1):

        if outputFileTemplate:
            fileName = outputFileTemplate % (alignmentNumber, str(l))
        else:
            fileName = listdir + str(alignmentNumber) + "." + str(l) + ".txt"

        writeList(alignmentNumber, l, hdf, fileName, chunkdir)

# def writeAllLists(alignmentNumber, hdf, chunkdir, listdir, outputFileTemplate=''):
#     if not listdir.endswith("/"):
#         listdir = listdir + "/"
#     for l in xrange(1, numberOfLists(alignmentNumber, hdf) +1):
#         writeList(alignmentNumber, l, hdf, listdir + str(alignmentNumber) + "." + str(l) + ".txt", chunkdir)
