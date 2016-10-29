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
    dst = hdf.create_table(root.segments, name, Segments)

#     state = None
#     buff = []
#     lastposSpecies = -1
#     lastposAln = -1

    prevState = None
    buff = []

    #for r in table.iterrows():
    for r in table.itersorted(table.cols.speciesPosition):
        (_, _, _, _, maxstate, maxP, chunk, alignmentPosition, alignmentNumber, speciesPosition) = r.fetch_all_fields()

        if speciesPosition == -1:
            continue # coord does not map to a speciesPosition in that species

        # first pos
        if prevState is None:
            segmentSpeciesStart = speciesPosition
            segmentAlnStart = alignmentPosition
            segmentSpeciesEnd = speciesPosition + 1
            segmentAlnEnd = alignmentPosition + 1
            prevState = maxstate
            continue
        
        # gap or change of max state (not evaluated in first loop)
        if speciesPosition > segmentSpeciesEnd or prevState != maxstate:
            buff.append((segmentAlnStart, segmentAlnEnd, segmentSpeciesStart, segmentSpeciesEnd, prevState))
            segmentSpeciesStart = speciesPosition
            segmentAlnStart = alignmentPosition

        # increment segment ends
        segmentSpeciesEnd = speciesPosition + 1
        segmentAlnEnd = alignmentPosition + 1        
        prevState = maxstate

        # write to table so we do not thrash memory
        if len(buff) > 10000:
            dst.append(buff)
            buff = []

    # add last segment and flush buffer
    buff.append((segmentAlnStart, segmentAlnEnd, segmentSpeciesStart, segmentSpeciesEnd, prevState))
    dst.append(buff) # write to table
        
        
#         if speciesPosition == -1:
#             continue # coord does not map to a speciesPosition in that species
# 
#         if state is None:
#             state = maxstate
#             startSpecies = speciesPosition
#             startAln = alignmentPosition
#             continue
# 
#         assert lastposSpecies < speciesPosition
# 
#         # start new segment if there is a gap (unless this is the first loop):
#         if speciesPosition > lastposSpecies + 1 and not lastposSpecies == -1:
#             buf.append((startAln, lastposAln+1, startSpecies, lastposSpecies+1, state))
#             startSpecies = speciesPosition
#             startAln = speciesAln
# 
#         lastposSpecies = speciesPosition
#         lastposAln = speciesAln
# 
#         # start new segment if there is a change to a new state:
#         if state != maxstate:                    
#             buf.append((startAln, lastposAln+1, startSpecies, lastposSpecies+1, state))
#             startSpecies = speciesPosition
#             startAln = speciesAln
#             state = maxstate
# 
#         if len(buff) > 10000: #write to table so we do not thrash memory
#             dst.append(buff)
#             buff = []
# 
#     buf.append((startAln, lastposAln+1, startSpecies, lastposSpecies+1, state))
#     dst.append(buff) # write to table

                    
