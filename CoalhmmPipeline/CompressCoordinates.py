from tables import *
from PipelineUtils import *

from time import time

def compressCoordinates(coords):
    chromosome = -1
    strand = 0
    alignmentNumber = -1
    start_sp = -1
    start_al = -1
    end_sp = -1
    end_al = -1
    count = 0
    total = 0
    
    # create tmp table in same group as coordinates
    if "_compress_tmp_" in coords._v_parent:
        print "error: temporary table allready exists"
        return
        
    # create table
    table = coords._v_file.createTable(coords._v_parent, "_compress_tmp_", SpeciesIntervalCoordinates, filters=Filters(complevel=0, complib='blosc', shuffle=True, fletcher32=False))
    buff = []


    #### this construct is about four times more efficient 
    chunk = 1000000
    data = []
    tabIdx = 0
    t = time()
    while True:
        data = coords.read(tabIdx, tabIdx+chunk)
        if not len(data):
            break
        print tabIdx, time() - t
        t = time()
        for a in data:
            x_chr, x_strand, x_sp, x_alno, x_al = a.tolist()
            tabIdx += 1
    #####################
            # is current pos NOT part of consequtive sequence of coordinates
            if chromosome != x_chr or \
               strand != x_strand or \
               alignmentNumber != x_alno or \
               abs(end_sp - x_sp) != 1 or \
               abs(end_al - x_al) != 1 or \
               (end_sp - start_sp)**2 > (x_sp - start_sp)**2 or \
               (end_al - start_al)**2 > (x_al - start_al)**2:
                #emit row
                if strand != 0: #test for first iteration
                    # output to buffer
                    buff.append((chromosome, strand, alignmentNumber, start_sp, end_sp, start_al, end_al))
                    if len(buff) >= 10:
                        # flush buffer
                        count = count + len(buff)
                        table.append(buff)
                        buff = []
#                        print "processed", count, "total length", total, "reduction", float(count)/float(total)

                # set new start parameters
                chromosome = x_chr
                strand = x_strand
                alignmentNumber = x_alno
                start_sp = x_sp
                start_al = x_al

            # update previous versions of parameters
            end_sp = x_sp
            end_al = x_al
            total = total +1

# 
#     for x in coords.iterrows():
#         (x_chr, x_strand, x_sp, x_alno, x_al) = x.fetch_all_fields()
# 
#         if chromosome != x_chr or \
#            strand != x_strand or \
#            alignmentNumber != x_alno or \
#            abs(end_sp - x_sp) != 1 or \
#            abs(end_al - x_al) != 1 or \
#            (end_sp - start_sp)**2 > (x_sp - start_sp)**2 or \
#            (end_al - start_al)**2 > (x_al - start_al)**2:
#             #emit row
#             if strand != 0: #test for first iteration
#                 buff.append((chromosome, strand, alignmentNumber, start_sp, end_sp, start_al, end_al))
#                 #print (chromosome, strand, alignmentNumber, start_sp, end_sp, start_al, end_al)
#                 if len(buff) >= 10:
#                     count = count + len(buff)
#                     table.append(buff)
#                     buff = []
# #                    print "processed", count, "total length", total, "reduction", float(count)/float(total)
#             
#             chromosome = x_chr
#             strand = x_strand
#             alignmentNumber = x_alno
#             start_sp = x_sp
#             start_al = x_al
#            
#         end_sp = x_sp
#         end_al = x_al
#         total = total +1
        
    # append last interval to buffer
    buff.append((chromosome, strand, alignmentNumber, start_sp, end_sp, start_al, end_al)) #assume number of rows in coords > 0
    print len(table)
    count = count + len(buff)
    # flush buffer
    if len(buff) > 0:
        table.append(buff)
        #print "processed", count, "total length", total, "reduction", float(count)/float(total)
#     print 'indexing'

    # remove old table and rename tmp table
    name = coords._v_name
    coords._f_remove()
    table._f_rename(name)
#     table.cols.alignmentNumber.createIndex()
#     table.cols.SpeciesPositionOnPlusStrandBegin.createIndex()
#     table.cols.SpeciesPositionOnPlusStrandEnd.createIndex()
#     table.cols.alignmentPositionBegin.createIndex()
#     table.cols.alignmentPositionEnd.createIndex()
#     print 'done'
            
#hdf = openFile("/tmp/julien.h5", "a")
#compressCoordinates(hdf.root.coordinates.human)













# import multiprocessing
# 
# def compressCoordinates(hdf, ingroup):    
# 
# 
# 
#     tabList = [hdf.root.coordinates[i] for i in ingroup]
# 
#     pool = multiprocessing.Pool(processes=4)
#     result = pool.apply_async(getCompressedData, tabList)
# 
#     for i in range(tabList):
#         coords, compressed = result.get()
# 
#         if "_compress_tmp_" in coords._v_parent:
#             print "error: temporary table allready exists"
#             return
# 
#         table = coords._v_file.createTable(coords._v_parent, "_compress_tmp_", SpeciesIntervalCoordinates,
#                                            filters=Filters(complevel=0, complib='blosc', shuffle=True, fletcher32=False))
#         buff = getCompressedData(coords)
#         table.append(buff)
# 
#         print 'indexing'
#         name = coords._v_name
#         coords._f_remove()
#         table._f_rename(name)
#         table.cols.alignmentNumber.createIndex()
#         table.cols.SpeciesPositionOnPlusStrandBegin.createIndex()
#         table.cols.SpeciesPositionOnPlusStrandEnd.createIndex()
#         table.cols.alignmentPositionBegin.createIndex()
#         table.cols.alignmentPositionEnd.createIndex()
#         print 'done'
# 
# 
# 
# def getCompressedData(coords):
#     chromosome = -1
#     strand = 0
#     alignmentNumber = -1
#     start_sp = -1
#     start_al = -1
#     end_sp = -1
#     end_al = -1
#     count = 0
#     total = 0
#     buff = []
# 
#     #### this construct is about four times more efficient 
#     chunk = 1000000
#     data = []
#     tabIdx = 0
#     t = time()
#     while True:
#         data = coords.read(tabIdx, tabIdx+chunk)
#         if not len(data):
#             break
#         print tabIdx, time() - t
#         t = time()
#         for a in data:
#             x_chr, x_strand, x_sp, x_alno, x_al = a.tolist()
#             tabIdx += 1
#     #####################
#             if chromosome != x_chr or \
#                strand != x_strand or \
#                alignmentNumber != x_alno or \
#                abs(end_sp - x_sp) != 1 or \
#                abs(end_al - x_al) != 1 or \
#                (end_sp - start_sp)**2 > (x_sp - start_sp)**2 or \
#                (end_al - start_al)**2 > (x_al - start_al)**2:
#                 #emit row
#                 if strand != 0: #test for first iteration
#                     buff.append((chromosome, strand, alignmentNumber, start_sp, end_sp, start_al, end_al))
# 
#                 chromosome = x_chr
#                 strand = x_strand
#                 alignmentNumber = x_alno
#                 start_sp = x_sp
#                 start_al = x_al
# 
#             end_sp = x_sp
#             end_al = x_al
#             total = total +1
#     
#     buff.append((chromosome, strand, alignmentNumber, start_sp, end_sp, start_al, end_al)) #assume number of rows in coords > 0
#     count = count + len(buff)
# 
#     return coords, buff
# 
# 
# 
#             
# #hdf = openFile("/tmp/julien.h5", "a")
# #compressCoordinates(hdf.root.coordinates.human)
#     
