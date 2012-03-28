from tables import *
import tables as pytables
from PipelineUtils import *
import numpy

from time import time


def mergeHDFfiles(sourceFileList, destFileName):

    newHDF = openFile(destFileName, "w")
    sourceHDFs = [openFile(s, 'r') for s in sourceFileList]
    merge([hdf.root for hdf in sourceHDFs], newHDF.root)
    for s in sourceHDFs:
        s.close()
    newHDF.close()


def merge(srcs, dst):

    assert dst._v_nchildren == 0
    for src in srcs:
#        print "processing", src._v_file.filename
        groups = []
        tables = []
        ##create groups and accumulate tables
        for g in src:
            if type(g) == pytables.group.Group:
                groups.append((g, dst))
            else:
                tables.append((g, dst))
                
        while len(groups) > 0:
            (gs, gd) = groups[0]
            groups = groups[1:]
            if gs._v_name not in gd:
                dst._v_file.createGroup(gd, gs._v_name)
                
            gdc = gd[gs._v_name]
            for g in gs:
                if type(g) == pytables.group.Group:
                    groups.append((g, gdc))
                else:
                    tables.append((g, gdc))

        #copy table data                    
        for (tab, grp) in tables:
#            print "copying", tab
            desttab = None
            if tab._v_name not in grp:
                desttab = dst._v_file.createTable(grp, tab._v_name, tab.description, filters=tab.filters)
            else:
                desttab = grp[tab._v_name]

            c = 1e7
            start, stop = 0, c
            size = len(tab)
            while start < size:
                desttab.append(tab.read(start, stop))
                start, stop = stop, stop + c

#             buff = []
#             t = time()
#             for x in tab:
#                 buff.append(x.fetch_all_fields())
# 
#                 if len(buff) >= 100000:
#                     print 'collect', time() - t
#                     t = time()
#                     desttab.append(numpy.array(buff, dtype=tab.description._v_dtype))
#                     print 'insert', time() - t
#                     t = time()
#                     buff = []
# 
#             if len(buff) > 0:
#                 desttab.append(numpy.array(buff, dtype=tab.description._v_dtype))
#                 buff = []

    for table in dst.coordinates:
        table.cols.alignmentNumber.createIndex()
        table.cols.SpeciesPositionOnPlusStrandBegin.createIndex()
        table.cols.SpeciesPositionOnPlusStrandEnd.createIndex()
        table.cols.alignmentPositionBegin.createIndex()
        table.cols.alignmentPositionEnd.createIndex()



#         #copy table data                    
#         for (tab, grp) in tables:
#             print "copying", tab
#             desttab = None
#             if tab._v_name not in grp:
#                 desttab = dst._v_file.createTable(grp, tab._v_name, tab.description, filters=tab.filters)
#             else:
#                 desttab = grp[tab._v_name]
#             buff = []
#             t = time()
#             for x in tab:
#                 buff.append(x.fetch_all_fields())
#                 
#                 if len(buff) >= 100000:
#                     print 'collect', time() - t
#                     t = time()
#                     desttab.append(numpy.array(buff, dtype=tab.description._v_dtype))
#                     print 'insert', time() - t
#                     t = time()
#                     buff = []
#                     
#             if len(buff) > 0:
#                 desttab.append(numpy.array(buff, dtype=tab.description._v_dtype))
#                 buff = []
                
#hdfs = openFile("/tmp/kasperTest1.1.h5", "r")
#hdfd = openFile("/tmp/test.h5", "w")

#merge([hdfs.root, hdfs.root], hdfd.root)
                
            
