from tables import *
from PipelineUtils import *

def createSegments(table):
    #check a few requirements
    #location should be /posteriors/$table
    if table._v_depth != 2:
        print "ERROR: table should be located at /posteriors/$table"
        return
    if table._v_parent._v_name != 'posteriors':
        print "ERROR: table should be located at /posteriors/$table"
        return
        
    root = table._v_parent._v_parent
    hdf = table._v_file
    name = table._v_name
    
    #test if group segments exists
    if "segments" not in root:
        hdf.createGroup(root, "segments")
        
    #test that a table with the same name does not allready exists
    if name in root.segments:
        print "ERROR: table named", name, "allready exists in /segments"
        return
        
    #create table
    dst = hdf.createTable(root.segments, name, Segments)

    processed = 0
    albegin = 0
    alend = 0
    spbegin = -123
    spend = spbegin
    state = -1
    buff = []
    for r in table.iterrows():
        (v0, v1, v2, v3, maxstate, maxP, chunk, alignmentPosition, alignmentNumber, speciesPosition) = r.fetch_all_fields()
        
#         if processed % 10000 == 0:
#             print "processed", processed
            
        processed = processed +1
        
        if state != maxstate or abs(speciesPosition - spend) != 1: #assumes case 1 2 3 2 case does not exist
            if state != -1: # simple test if we are in a valid interval (not in first run)
                buff.append((albegin, alend+1, spbegin, spend+1, state)) #emit interval (+1 to make interval [,[ instead of [,])
                
            spbegin = speciesPosition
            albegin = alignmentPosition
            state = maxstate
        
        spend = speciesPosition
        alend = alignmentPosition
        
        if len(buff) > 10000: #write to table so we do not thrash memory
            dst.append(buff)
            buff = []
            
    buff.append((albegin, alend+1, spbegin, spend+1, state)) #add last interval
    dst.append(buff) #write to table    
